#!/bin/env perl
########################################################################
#
#Author: Haisi Yi
#Contact: yihaisi@gmail.com
#
########################################################################

use strict;
#use warnings;
#use Memory::Stats;
#use Devel::Size qw(size total_size);
use Getopt::Long;
use List::Util qw(sum product max min);
#use String::Approx qw(amatch adist);
use Data::Dumper;

########################################################################


#####################
#### MAIN begins ####
#####################

#my $stats = Memory::Stats->new; # monitor the memory usage
#$stats->start;

print "Command line: $0 @ARGV\n"; # print the command line

#### set default values of parameters
my @OPT_input          = ();
my @OPT_idxInput       = ();
my $OPT_output         = '';
my $OPT_idxLocation    = '';
my $OPT_qType          = 33;

my $OPT_adapterSeq1    = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC';
my $OPT_adapterSeq2    = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT';
my $OPT_useAdapter     = 8;
my $OPT_idx2direction  = 'ff'; 
my $OPT_minTrainQual   = 5;

my $OPT_priorP         = 'auto';
my $OPT_naiveUse       = 5;
my $OPT_useOVA         = 150;
my $OPT_confidence     = 0.95;
my $OPT_throwLine      = 10;

my $MLminQuality   = 10;
my $ML             = 0;
my $MLZ1           = 0.01;

my $help           = 0;
my $dev            = 0;

my $usage = "
About:
This program is for demultiplexing the multiplexed Illumina sequencing data.
Contact: Haisi Yi <yihaisi\@ihb.ac.cn>

Input/Ouput Options:
-i    the fastq file(s) of multiplexed Illumina reads
-j    the fastq file(s) of index reads
-o    the output directory in which the demultiplexed fastq files will be put
-x    the file of sample-index list
-q    the quality score type(phred33 or phred64) [$OPT_qType]

Taining Data Extraction Options
-a    the pre-index1 adapter sequence [$OPT_adapterSeq1]
-b    the pre-index2 adapter sequence [$OPT_adapterSeq2]
-u    use last N bases of the pre-index sequence in the search [$OPT_useAdapter]
-d    the relative direction of the index 2 and its upstream adapter sequence (ff or fr) [$OPT_idx2direction]
-n    the minimum quality score of the index bases to accept in the adapter searching [$OPT_minTrainQual]

Bayes Options
-p    turn on/off the auto inference of prior probability(auto/infer) [$OPT_priorP]
-f    the minimum number of evidences of a feature to be used in Naive Bayes [$OPT_naiveUse] 
-v    the maximum occurrence of a barcode to use the one-versus-all-but-one technique [$OPT_useOVA]
-c    the minimum P to be trusted [$OPT_confidence]
-l    the minimum occurence frequency of a sequenced index be considered in the Bayes Module [$OPT_throwLine]

--help    this help information
--dev    for develop use only
";
die $usage if @ARGV == 0;
####

#### get parameters
GetOptions('i=s{1,2}'   => \@OPT_input,
	       'j=s{1,2}'   => \@OPT_idxInput,
		   'o=s'   => \$OPT_output,
		   'x=s'   => \$OPT_idxLocation,
		   'a=s'   => \$OPT_adapterSeq1,
		   'b=s'   => \$OPT_adapterSeq2,
		   'd=s'   => \$OPT_idx2direction,
		   'u=i'   => \$OPT_useAdapter,
		   'n=i'   => \$OPT_minTrainQual,
		   'l=f'   => \$OPT_throwLine,
		   'p=s'   => \$OPT_priorP,
		   'f=i'   => \$OPT_naiveUse,
		   'v=i'   => \$OPT_useOVA,
		   'c=f'   => \$OPT_confidence,
		   'ML!'   => \$ML,
		   't=f'   => \$MLminQuality,
		   'F=f'   => \$MLZ1,
		   'q=i'   => \$OPT_qType,
		   'help!' => \$help,
           'dev!'  => \$dev);
####
		   
#$stats->checkpoint("before data process.");


#### check all the parameters 
my $indexTYPE = 0;
&checkParameters();
####


#### parse the index sheet file
my %hCodeToName = ();
my %hPctToCode = ();
my %hCodeToClass = ();
my %outhandle = ();
my %outIdxHandle = ();
my %hIndex1WithAdapter = ();
my %hIndex2WithAdapter = ();
my $inferP = 0;

if($indexTYPE == 2){
	&parseIndexSheet();
}
elsif($indexTYPE == 1){
	&parseIndexSheet_sin();
}
####

my %hClassToCode = reverse%hCodeToClass;
my @aIndexDual = keys(%hCodeToClass);
my $indexTotal = scalar @aIndexDual;
my $indexOneLength = 0;
if($indexTYPE == 2){
	$indexOneLength = (length((keys%hCodeToClass)[0]))/2;
}
elsif($indexTYPE == 1){
	$indexOneLength = length((keys%hCodeToClass)[0]);
}

#$stats->checkpoint("after index read");


#### seatch the train data set
my %hBothSum = ();
my @Train = ();
my %hAllBarcode = ();
my %hDirectBase = ();
my $sumEvidenceBoth = 0;

if($indexTYPE == 2){
	&searchTrainData();
}
elsif($indexTYPE == 1){
	&searchTrainData_sin();
}
####

#$stats->checkpoint("after search evidence");

my ($levelCount, $barcodeOneLength) = (0, 0);
if($indexTYPE == 2){
	$barcodeOneLength = sqrt($#Train+0.25)-0.5;
	$levelCount = $barcodeOneLength+1;
}
elsif($indexTYPE == 1){
	$barcodeOneLength = sqrt(2*($#Train+1)+0.25)-0.5;
	$levelCount = $barcodeOneLength;
}
my $barcodeLength = 2*$barcodeOneLength;

#### infer the prior probability
&inferPriorP();
####


#### smoothing and P calculation
my @TrainSec = ();
my @TrainSecOVA = ();

&computeTrainSec();
&computeTrainSecOVA();
####

#$stats->checkpoint("after P calc");

#### calculate the Ps for every barcode
my @FeatureForSubstr = &buildFeatureForSbustr();
my @featurePostion = &makeFeaturePostion();

my %hForDecision = ();
&computePostP();
####

#$stats->checkpoint("after secondary calc");

#print total_size(\@Train),"\t",total_size(\@TrainSec),"\t",total_size(\%TrainSecOVA),"\n"; # for test only

#### dump the big data structures (for develop use only)
&dumpBigData() if $dev == 1;
####

#### clear big data structure for memory
undef @Train;
undef @TrainSec;
undef @TrainSecOVA;
####

#### decision by threshold OPT_throwLine and OPT_confidence
my %hDecisionBase = ();
my %summary = ();
my %summaryTotal = ();
&makeDicisionAndSTATS();
####

#$stats->checkpoint("after decision base built");

#### clear big data structure for memory
undef %hForDecision;
undef %hAllBarcode;
####

#### dump the assignment decisions (for develop use only)
if($dev == 1){
	open DECI,'>>',"$OPT_output/decision.dump";
	print DECI Dumper(%hDecisionBase); 
	close DECI;
}
####

#### start the file reading and demultiplexing
if($indexTYPE == 2){
	&outputDATA();
}
elsif($indexTYPE == 1){
	&outputDATA_sin();
}
####

#### output the summary file
&outputSummary();
####

#$stats->checkpoint("at last");

print &checkTime()."Mission complete!\n";

#$stats->stop;
#$stats->report;

exit(0);

###################
#### MAIN ends ####
###################
		   
		   
##########################################################
		   
		   
#########################
#### FUNCTIONS begin ####
#########################

#### tools ####

#check parameters
sub checkParameters {
	die $usage if $help==1;
	map{s/\/$//}(@OPT_input,@OPT_idxInput,$OPT_output);

	die "Error: Only one or two fastq files are acceptable for the input after -i.\n" if @OPT_input != 2 and @OPT_input != 1;
	die "Error: Only one or two fastq files are acceptable for the input after -j.\n" if @OPT_idxInput != 2 and @OPT_idxInput !=1;
	$indexTYPE = scalar @OPT_idxInput; #### detect single-index or dual-index
	die "Error: No index file after -x !\n" if ! $OPT_idxLocation;
	die "Error: No output directory after -o !\n" if ! $OPT_output;
	die "Error: Illegal value '$OPT_priorP' for -p.\n" if $OPT_priorP ne 'auto' and $OPT_priorP ne 'infer';
	die "Error: Illegal value '$OPT_idx2direction' for -d.\n" if $OPT_idx2direction ne 'ff' and $OPT_idx2direction ne 'fr';
		
	for my $file (@OPT_input, @OPT_idxInput, $OPT_idxLocation){
		die "Error: $file is not an available text file. $!.\n" if ! -T $file;
	}
	if(!-e $OPT_output){
		mkdir("$OPT_output",0755) or die "Error: Can not make directory $OPT_output. $!.\n";
	}
	else{
		die "Error: $OPT_output already exists. $!.\n";
	}
}

sub reverseComplement {
	my $a = shift;
	$a =~ tr{AaTtCcGgNn}{TtAaGgCcNn};
	return scalar reverse $a;
}

sub makeKmerArray {
	my ($string1, $string2) = (shift, shift);
	my $l = length$string1;
	my @kArray = ($string1.$string2);
	foreach my $k (reverse 1 .. $l){
		push(@kArray, substr($string1, $_, $k), substr($string2, $_, $k)) foreach ( 0 .. ($l - $k));
	}
	return @kArray;
}

sub makeKmerArray_sin {
	my $string = shift;
	my $l = length$string;
	my @kArray = ();
	foreach my $k (reverse 1 .. $l){
		push(@kArray, substr($string, $_, $k)) foreach (0 .. ($l - $k));
	}
	return @kArray;
}

sub calcML { # quality, barcode, index
	my $q = shift;
	my $b = shift;
	my $i = shift;
	my $P = 1;
	foreach my $n (0 .. (length$i)-1){
		#my $c = ord(substr($q,$n,1)) - $OPT_qType;	
		#next if $c < 3;
		$P *= abs(substr($b,$n,1) eq substr($i,$n,1)) - 10**(-(ord(substr($q,$n,1))-$OPT_qType)/10);
	}
	return $P;
}

sub checkTime {
	my $t = times;
	my $m = int $t/60;
	my $s = $t%60;
	return "${m}m${s}s\t";
}

##################################################################

#### process ####

#### parse the index file and prepare the output file handles
sub parseIndexSheet {
	print &checkTime()."Reference Index sheet parse begins.\n";
	
	open INDEX,"<","$OPT_idxLocation" or die "Can not read index file $OPT_idxLocation: $!\n";
	my $id = 0;
	
	while(my $indline = <INDEX>){
		next if $indline =~ m/^\s*$/ or $indline =~ m/^#/;
		my @index = split(/\s+/, $indline);
		die("Incorrect index sequence in line $. of $ARGV.\n") if $index[1]!~/^[ATCGNatcgn]+$/ or $index[2]!~/^[ATCGNatcgn]+$/;
		
		$hCodeToName{$index[1].$index[2]} = $index[0];
		$hPctToCode{$index[1].$index[2]} = $index[3];
		$hCodeToClass{$index[1].$index[2]} = ++$id;
		$inferP = 1 unless $index[3] > 0;
		$hIndex1WithAdapter{substr($OPT_adapterSeq1,-$OPT_useAdapter,$OPT_useAdapter).$index[1]} = $index[1];
		$hIndex2WithAdapter{substr($OPT_adapterSeq2,-$OPT_useAdapter,$OPT_useAdapter).(($OPT_idx2direction eq "fr") ? &reverseComplement($index[2]) : $index[2])} = $index[2];
	
		!-e "$OPT_output/$index[0].read1.fastq"
			? open($outhandle{$id."R1"},">>","$OPT_output/$index[0].read1.fastq") : die $!;
		!-e "$OPT_output/$index[0].read1.untrusted.fastq"
			? open($outhandle{$id."UNR1"},">>","$OPT_output/$index[0].read1.untrusted.fastq") : die $!;
		!-e "$OPT_output/$index[0].read2.fastq"
			? open($outhandle{$id."R2"},">>","$OPT_output/$index[0].read2.fastq") : die $!;
		!-e "$OPT_output/$index[0].read2.untrusted.fastq"
			? open($outhandle{$id."UNR2"},">>","$OPT_output/$index[0].read2.untrusted.fastq") : die $!;
		!-e "$OPT_output/$index[0].index1.fastq"
			? open($outIdxHandle{$id."I1"},">>","$OPT_output/$index[0].index1.fastq") : die $!;
		!-e "$OPT_output/$index[0].index1.untrusted.fastq"
			? open($outIdxHandle{$id."UNI1"},">>","$OPT_output/$index[0].index1.untrusted.fastq") : die $!;
		!-e "$OPT_output/$index[0].index2.fastq"
			? open($outIdxHandle{$id."I2"},">>","$OPT_output/$index[0].index2.fastq") : die $!;
		!-e "$OPT_output/$index[0].index2.untrusted.fastq"
			? open($outIdxHandle{$id."UNI2"},">>","$OPT_output/$index[0].index2.untrusted.fastq") : die $!;
	}
	close INDEX;

	print &checkTime()."Reference index sheet parse finished and output files prepared.\n";
}
####

sub parseIndexSheet_sin {
	print &checkTime()."Reference Index sheet parse begins.\n";

	open INDEX,"<","$OPT_idxLocation" or die "Can not read index file $OPT_idxLocation: $!\n";
	my $id = 0;

	while(my $indline = <INDEX>){
		next if $indline =~ m/^\s*$/ or $indline =~ m/^#/;
		my @index = split(/\s+/, $indline);
		die("Incorrect index sequence in line $. of $ARGV.\n") if $index[1]!~/^[ATCGNatcgn]+$/;
		
		$hCodeToName{$index[1]} = $index[0];
		$hPctToCode{$index[1]} = $index[2];
		$hCodeToClass{$index[1]} = ++$id;
		$inferP = 1 unless $index[2] > 0;
		$hIndex1WithAdapter{substr($OPT_adapterSeq1,-$OPT_useAdapter,$OPT_useAdapter).$index[1]} = $index[1];
	
		foreach my $n (1 .. scalar@OPT_input){
			!-e "$OPT_output/$index[0].read$n.fastq"
				? open($outhandle{$id."R".$n},">>","$OPT_output/$index[0].read$n.fastq") : die $!;
			!-e "$OPT_output/$index[0].read$n.untrusted.fastq"
				? open($outhandle{$id."UNR".$n},">>","$OPT_output/$index[0].read$n.untrusted.fastq") : die $!;
		}
		!-e "$OPT_output/$index[0].index.fastq"
			? open($outIdxHandle{$id."I"},">>","$OPT_output/$index[0].index.fastq") : die $!;
		!-e "$OPT_output/$index[0].index.untrusted.fastq"
			? open($outIdxHandle{$id."UNI"},">>","$OPT_output/$index[0].index.untrusted.fastq") : die $!;
	}
	close INDEX;

	print &checkTime()."Reference index sheet parse finished and output files prepared.\n";
}

#### analysis the contaminant sequences containing barcode bases
sub searchTrainData {
	print &checkTime()."Train Set search begins. (This process may take a while, please wait.)\n";
	
	my @inHandle = ();
	open($inHandle[$_], '<', "$OPT_input[$_]") or die $! foreach 0 .. $#OPT_input;
	my @idxHandle = ();
	open($idxHandle[$_], '<', "$OPT_idxInput[$_]") or die $! foreach 0 .. $#OPT_idxInput;

	my $minChar = $OPT_qType+$OPT_minTrainQual;
	my $counter = 0;
	SEARCH: while(readline(*{$inHandle[0]})){
		my $read1 = readline(*{$inHandle[0]});
		readline(*{$inHandle[0]});
		my $qual1 = readline(*{$inHandle[0]});

		readline(*{$inHandle[1]});
		my $read2 = readline(*{$inHandle[1]});
		readline(*{$inHandle[1]});
		my $qual2 = readline(*{$inHandle[1]});

		readline(*{$idxHandle[0]});
		my $barcode1 = readline(*{$idxHandle[0]});
		readline(*{$idxHandle[0]});
		readline(*{$idxHandle[0]});

		readline(*{$idxHandle[1]});
		my $barcode2 = readline(*{$idxHandle[1]});
		readline(*{$idxHandle[1]});
		readline(*{$idxHandle[1]});

		$counter++; # count for direct base

		chomp ($barcode1,$barcode2);
		$hAllBarcode{$barcode1.$barcode2} ++; # build hash of all existing Barcodes

		#### search for evidences
		my ($find1, $find2, $pos1, $pos2) = (0, 0, 0, 0);
		foreach my $searchSeq (keys%hIndex1WithAdapter){
			if($pos1 = index($read1, $searchSeq) + 1){ # find adapters contain barcodes
				$find1 = $hIndex1WithAdapter{$searchSeq};
				last;
			}
		}
		next SEARCH unless $find1;
		foreach my $searchSeq (keys%hIndex2WithAdapter){
			if($pos2 = index($read2, $searchSeq) + 1){
				$find2 = $hIndex2WithAdapter{$searchSeq};
				last;
			}
		}
		next SEARCH unless $find2;
		#### search end

		#### build module
		my $find = $find1.$find2;
		if(exists$hCodeToName{$find}){
			foreach ($pos1-1+$OPT_useAdapter .. $pos1+$OPT_useAdapter+$indexOneLength){
				next SEARCH if vec($qual1, $_, 8) < $minChar;
			}
			foreach ($pos2-1+$OPT_useAdapter .. $pos2+$OPT_useAdapter+$indexOneLength){
				next SEARCH if vec($qual2, $_, 8) < $minChar;
			}
			$hDirectBase{$counter} = $hCodeToClass{$find};
			$hBothSum{$hCodeToClass{$find}} ++;
			my $tn = 0;
			$Train[$tn++]{$_}{$hCodeToClass{$find}} ++ foreach &makeKmerArray($barcode1,$barcode2); # Naive Bayes Module
		}
		#### build end
	}
	close $_ for @inHandle,@idxHandle;
	
	print &checkTime()."Train Set search finished.\n";
	$sumEvidenceBoth = sum(values%hBothSum);
	print "Total evidence: $sumEvidenceBoth\n";
	print "\t$hClassToCode{$_}: $hBothSum{$_}\n" for keys%hBothSum;
}

sub searchTrainData_sin {
	print &checkTime()."Train Set search begins. (This process may take a while, please wait.)\n";

	my @inHandle = ();
	open($inHandle[0], '<', "$OPT_input[0]") or die $!;
	my @idxHandle = ();
	open($idxHandle[0], '<', "$OPT_idxInput[0]") or die $!;

	my $minChar = $OPT_qType+$OPT_minTrainQual;
	my $counter = 0;
	SEARCH: while(readline(*{$inHandle[0]})){
		my $read = readline(*{$inHandle[0]});
		readline(*{$inHandle[0]});
		my $qual = readline(*{$inHandle[0]});

		readline(*{$idxHandle[0]});
		my $barcode = readline(*{$idxHandle[0]});
		readline(*{$idxHandle[0]});
		readline(*{$idxHandle[0]});

		$counter++; # count for direct base

		chomp($barcode);
		$hAllBarcode{$barcode} ++; # build hash of all existing Barcodes

		#### search for evidences
		my ($find, $pos) = (0, 0);
		foreach my $searchSeq (keys%hIndex1WithAdapter){
			if($pos = index($read, $searchSeq) + 1){ # find adapters contain barcodes
				$find = $hIndex1WithAdapter{$searchSeq};
				last;
			}
		}
		next SEARCH unless $find;
		#### search end

		#### build module
		foreach ($pos - 1 + $OPT_useAdapter .. $pos + $OPT_useAdapter + $indexOneLength){
			next SEARCH if vec($qual, $_, 8) < $minChar;
		}
		$hDirectBase{$counter} = $hCodeToClass{$find};
		$hBothSum{$hCodeToClass{$find}}++;
		my $tn = 0;
		$Train[$tn++]{$_}{$hCodeToClass{$find}}++ foreach &makeKmerArray_sin($barcode); # Naive Bayes Module
		#### build end
	}
	close $_ for @inHandle,@idxHandle;
	
	print &checkTime()."Train Set search finished.\n";
	$sumEvidenceBoth = sum(values%hBothSum);
	print "Total evidence: $sumEvidenceBoth\n";
	print "\t$hClassToCode{$_}: $hBothSum{$_}\n" for keys%hBothSum;
}

####

#### prior Probability inference
sub inferPriorP {
	if($inferP == 1 or $OPT_priorP eq "infer"){
		foreach my $idx (@aIndexDual){
			$hPctToCode{$idx} = 1;
			foreach my $k (keys %hAllBarcode) {
				if($indexTYPE == 2){
					$hPctToCode{$idx}+=$hAllBarcode{$k} 
					if substr($idx,0,$indexOneLength) eq substr($k,0,$indexOneLength) 
					and substr($idx,$indexOneLength,$indexOneLength) eq substr($k,$barcodeOneLength,$indexOneLength);
				}
				elsif($indexTYPE == 1){
					$hPctToCode{$idx}+=$hAllBarcode{$k}
					if substr($idx,0,$indexOneLength) eq substr($k,0,$indexOneLength);
				}
			}
			print "Warning: prior probability inference on index '$idx' may be unreliable!\n" if $hPctToCode{$idx} < 5;
		}
	}
	my $sumPorpotion = sum(values%hPctToCode);
	$hPctToCode{$_} /= $sumPorpotion foreach keys%hPctToCode;

	print &checkTime()."Prior probability inference complete.\n";
}
####

sub computeTrainSec {	
	print &checkTime()."Bayes module: computation of P of each feature begins.\n";
	foreach my $feature (0 .. $#Train) {
		foreach my $k (keys%{$Train[$feature]}) {
			foreach my $idx (values%hCodeToClass){
				$TrainSec[$feature]{$k}{$idx} =
					($Train[$feature]{$k}{$idx} + (($hBothSum{$idx}+1/$indexTotal)/($sumEvidenceBoth+1))/$indexTotal )
					/
					($hBothSum{$idx} + (scalar keys%{$Train[$feature]})*((($hBothSum{$idx}+1/$indexTotal)/($sumEvidenceBoth+1))/$indexTotal) );
			}
		}
	}
	print &checkTime()."Bayes module: computation of P of each feature finished.\n";
}

sub computeTrainSecOVA {
	print &checkTime()."Bayes module: computation of OVA P of each feature begins.\n";
	foreach my $feature (0 .. $#Train) {
		foreach my $k (keys%{$Train[$feature]}) {
			foreach my $idx (values%hCodeToClass){
				$TrainSecOVA[$feature]{$k}{$idx} =
					((sum(values%hBothSum)-$hBothSum{$idx} + (scalar keys%{$Train[$feature]})*((($hBothSum{$idx}+1/$indexTotal)/($sumEvidenceBoth+1))/$indexTotal) )
					/
					(sum(values%{$Train[$feature]{$k}})-$Train[$feature]{$k}{$idx} + (($hBothSum{$idx}+1/$indexTotal)/($sumEvidenceBoth+1))/$indexTotal ))
					*
					(($Train[$feature]{$k}{$idx} + (($hBothSum{$idx}+1/$indexTotal)/($sumEvidenceBoth+1))/$indexTotal )
					/
					($hBothSum{$idx} + (scalar keys%{$Train[$feature]})*((($hBothSum{$idx}+1/$indexTotal)/($sumEvidenceBoth+1))/$indexTotal) ));
			}
		}
	}
	print &checkTime()."Bayes module: computation of OVA P of each feature finished.\n";
}

####

sub makeFeaturePostion {
	my @aPos = ();
	if($indexTYPE == 2){
		push(@aPos, [0, 0, 1]); # featureBegin featureEnd featureCount
		push(@aPos, [$_*($_-1)+1, $_*($_+1), 2*$_]) for (1 .. $barcodeOneLength);
	}
	elsif($indexTYPE == 1){
		push(@aPos, [0.5*($_-1)*$_, 0.5*$_*($_+1)-1, $_]) for (1 .. $barcodeOneLength);
	}
    return @aPos;
}

sub buildFeatureForSbustr {
	my $x = 0;
	my @FeatureForSubstr = ();
	if($indexTYPE == 2){
		push(@FeatureForSubstr, [0, $barcodeLength]);
		foreach my $L (0 .. $barcodeOneLength - 1) {
			foreach my $n (1 .. ($x += 2)) {
				push(@FeatureForSubstr, [int(($n - 1) / 2)+($n % 2 ? 0 : $barcodeOneLength), $barcodeOneLength-$L]);
			}
		}
	}
	elsif($indexTYPE == 1){
		foreach my $L (0 .. $barcodeOneLength - 1) {
			foreach my $n (0 .. $x++) {
				push(@FeatureForSubstr, [$n, $barcodeOneLength-$L]);
			}
		}
	}
	return @FeatureForSubstr;
}

sub computePostP {
	print &checkTime()."Bayes module: P calculation for every barcode begins.\n";
		
	foreach my $barcode (keys%hAllBarcode){
		my %P = ();
		my @useFeature = ();
		#my $L = 0; # for test only
		foreach my $l (0 .. $levelCount-1){
			#$L = $l; # for test only
			my @aNowUse = grep{
				my $K = substr($barcode, $FeatureForSubstr[$_][0], $FeatureForSubstr[$_][1]);
				exists($Train[$_]{$K}) and sum(values%{$Train[$_]{$K}}) >= $OPT_naiveUse;
			}($featurePostion[$l][0] .. $featurePostion[$l][1]);
			push(@useFeature, @aNowUse);
			last if scalar(@aNowUse) == $featurePostion[$l][2];
		}

		if($hAllBarcode{$barcode} > $OPT_useOVA){
			foreach my $idx (@aIndexDual) {
				$P{$idx} = 
					product(
						map({$TrainSec[$_]{substr($barcode, $FeatureForSubstr[$_][0], $FeatureForSubstr[$_][1])}{$hCodeToClass{$idx}}}@useFeature),
						$hPctToCode{$idx}
					);
			}
		}
		else{
			foreach my $idx (@aIndexDual) {
				$P{$idx} = 
					product(
						map({$TrainSecOVA[$_]{substr($barcode, $FeatureForSubstr[$_][0], $FeatureForSubstr[$_][1])}{$hCodeToClass{$idx}}}@useFeature),
						$hPctToCode{$idx}
					);
			}
		}
		$hForDecision{(sort{$P{$b} <=> $P{$a}}keys%P)[0]}{$barcode} = max(values%P) / (sum(values%P) || 1);
		#print $barcode,"\t",$L,"\t",(sort{$P{$b} <=> $P{$a}}keys%P)[0],"\t",(max(values%P)/sum(values%P)),"\n"; # for test only
		#print ("$barcode\t$L\t".($P{'GACGATTGACGGCG'}/sum(values%P))."\t".((sort{$P{$b} <=> $P{$a}}keys%P)[0])."\n") ; # for test only
	}
		
	print &checkTime()."Bayes module: P calculation for every barcode finished.\n";
}
####

sub dumpBigData {
	open Ev,'>>',"$OPT_output/evidence_for_bayesNaive.dump";
	print Ev Dumper(@Train);
	close Ev;
	
	open PP,">>","$OPT_output/priorP.dump";
	print PP Dumper(%hPctToCode);
	close PP;
	
	open Ba,'>>',"$OPT_output/bayesNaive.dump";
	print Ba Dumper(@TrainSec);
	close Ba;

	open DB,'>>',"$OPT_output/DirectDemultiplex.dump";
	print DB Dumper(%hDirectBase);
	close DB;

	open ForD,'>>',"$OPT_output/forDecision.dump";
	print ForD Dumper(%hForDecision);
	close ForD;
}

sub makeDicisionAndSTATS {
	print &checkTime()."Assignment of barcode-to-sample begins.\n";
	my $line = 0;
	if($OPT_throwLine >= 1 or $OPT_throwLine == 0){
		$line = $OPT_throwLine;
	}
	elsif($OPT_throwLine < 1 and $OPT_throwLine > 0){
		my $sumTH = $OPT_throwLine * sum(values(%hAllBarcode));
		my $level = 0;
		foreach my $c (sort{$b <=> $a} values%hAllBarcode){
			$line = $c;
			$level += $line;
			last if $level >= $sumTH;
		}
	}

	foreach my $idx (keys%hForDecision) {
		foreach my $barcode (keys%{$hForDecision{$idx}}){
			$summary{$idx}{'total'} += $hAllBarcode{$barcode};
			$summaryTotal{'total'} += $hAllBarcode{$barcode};
			if($hAllBarcode{$barcode} >= $line && $hForDecision{$idx}{$barcode} >= $OPT_confidence){
				$hDecisionBase{$barcode} = $hCodeToClass{$idx};
				$summary{$idx}{'trusted'} += $hAllBarcode{$barcode};
				$summaryTotal{'trusted'} += $hAllBarcode{$barcode};
			}
			else{
				$hDecisionBase{$barcode} = $hCodeToClass{$idx}.'UN';
				$summary{$idx}{'untrusted'} += $hAllBarcode{$barcode};
				$summaryTotal{'untrusted'} += $hAllBarcode{$barcode};
			}
		}
	}
	print &checkTime()."Assignment of barcode-to-sample finished.\n";
}
####

####
sub outputSummary {
	open(SUM,">>","$OPT_output/summary.txt") or die $!;
	print SUM "Sample\t".
			  "Total\t".
			  "Total%\t".
			  "Trusted\t".
			  "Trusted%\t".
			  "Untrusted\t".
			  "Untrusted%\n".
			  "----------------------------\n";
	foreach my $sample (sort{$summary{$b}{'total'} <=> $summary{$a}{'total'}}keys%summary){
		print SUM $hCodeToName{$sample}."\t".
				  $summary{$sample}{'total'}."\t".
				  (sprintf "%.3f%%",100*$summary{$sample}{total}/$summaryTotal{total}).
				  "\t".$summary{$sample}{'trusted'}."\t".
				  (sprintf "%.3f%%",100*$summary{$sample}{trusted}/$summary{$sample}{total}).
				  "\t".$summary{$sample}{'untrusted'}."\t".
				  (sprintf "%.3f%%",100*$summary{$sample}{untrusted}/$summary{$sample}{total}).
				  "\n";
	}
	print SUM "----------------------------\n".
			  "TOTOL NUMBER: $summaryTotal{'total'}\n".
			  "Trusted: $summaryTotal{'trusted'} ".(sprintf "%.3f%%",100*$summaryTotal{'trusted'}/$summaryTotal{'total'})."\n".
			  "Untrusted: $summaryTotal{'untrusted'} ".(sprintf "%.3f%%",100*$summaryTotal{'untrusted'}/$summaryTotal{'total'})."\n";
	close SUM;

	print &checkTime()."Summary file writed.\n";
}
####

sub outputDATA_sin {
	print &checkTime()."Output begins.\n";
	
	my @inHandle = ();
	open($inHandle[$_], '<', "$OPT_input[$_]") or die $! foreach 0 .. $#OPT_input;
	open(my $idxHandle, '<', "$OPT_idxInput[0]") or die $!;
	
	my $counter = 0;

	while(1){
		my (@R, @I) = ((), ());
		foreach my $n (0 .. $#OPT_input){
			push(@{$R[$n]}, scalar readline(*{$inHandle[$n]})) foreach 0 .. 3;
		}
		push(@I, scalar readline($idxHandle)) foreach 0 .. 3;
		
		$counter ++;
		
		my $barcode = $I[1];
		chomp $barcode;
		
		my $out = '';
		if(exists$hDirectBase{$counter}){
			$out = $hDirectBase{$counter};
			my ($n, $u) = $hDecisionBase{$barcode} =~ m/(\d+)(.*)/;
			my $ifTrust = $u eq 'UN' ? 'untrusted' : 'trusted';
			$summary{$hClassToCode{$n}}{$ifTrust}--;
			$summaryTotal{$ifTrust}--;
			$summary{$hClassToCode{$out}}{'trusted'}++;
			$summaryTotal{'trusted'}++;
		}else{
			$out = $hDecisionBase{$barcode};
		}
		foreach my $n (1 .. scalar@OPT_input){
			my $outRead = $outhandle{$out."R".$n};
			print $outRead @{$R[$n-1]};
		}
		my $outIdx = $outIdxHandle{$out."I"};
		print $outIdx @I;
		
		last if eof;
	}
	close $_ foreach @inHandle,$idxHandle;
	
	print &checkTime()."Output finished.\n";
}

sub outputDATA {

print &checkTime()."Output begins.\n";

#### open input file handles
my @inHandle = ();
open $inHandle[$_], '<', "$OPT_input[$_]" or die $! foreach 0 .. $#OPT_input;
my @idxHandle = ();
open $idxHandle[$_], '<', "$OPT_idxInput[$_]" or die $! foreach 0 .. $#OPT_idxInput;
####

#### read files and demultiplexing
my $counter = 0;

if($ML == 0){ # turn off ML algorithm
	while(1){
		my (@Fr1,@Fr2,@Ir1,@Ir2) = ((),(),(),());
		push(@Fr1, scalar readline(*{$inHandle[0]})) foreach 0 .. 3;
		push(@Fr2, scalar readline(*{$inHandle[1]})) foreach 0 .. 3;
		push(@Ir1, scalar readline(*{$idxHandle[0]})) foreach 0 .. 3;
		push(@Ir2, scalar readline(*{$idxHandle[1]})) foreach 0 .. 3;
	
		$counter++;
	
		my ($barcode1,$barcode2) = ($Ir1[1],$Ir2[1]);
		chomp ($barcode1, $barcode2);
	
		my $out = '';
		if(exists$hDirectBase{$counter}){
			$out = $hDirectBase{$counter};
			my ($n, $u) = $hDecisionBase{$barcode1.$barcode2} =~ m/(\d+)(.*)/;
			my $ifTrust = $u eq 'UN' ? 'untrusted' : 'trusted';
			$summary{$hClassToCode{$n}}{$ifTrust}--;
			$summaryTotal{$ifTrust}--;
			$summary{$hClassToCode{$out}}{'trusted'}++;
			$summaryTotal{'trusted'}++;
		}else{
			$out = $hDecisionBase{$barcode1.$barcode2};
		}

		my $outRead1 = $outhandle{$out."R1"};
		my $outRead2 = $outhandle{$out."R2"};
		my $outIdx1 = $outIdxHandle{$out."I1"};
		my $outIdx2 = $outIdxHandle{$out."I2"};
		print $outRead1 @Fr1;
		print $outRead2 @Fr2;
		print $outIdx1 @Ir1;
		print $outIdx2 @Ir2;
		
		last if eof;
	}
}
elsif($ML == 1){ # turn on ML algorithm
	while(1){
		my (@Fr1,@Fr2,@Ir1,@Ir2) = ((),(),(),());
		push(@Fr1, scalar readline(*{$inHandle[0]})) foreach 0 .. 3;
		push(@Fr2, scalar readline(*{$inHandle[1]})) foreach 0 .. 3;
		push(@Ir1, scalar readline(*{$idxHandle[0]})) foreach 0 .. 3;
		push(@Ir2, scalar readline(*{$idxHandle[1]})) foreach 0 .. 3;
	
		$counter++;
	
		my($barcode1,$barcode2,$quality1,$quality2) = ($Ir1[1],$Ir2[1],$Ir1[3],$Ir2[3]);
		chomp ($barcode1,$barcode2,$quality1,$quality2);
		my $barcode = $barcode1.$barcode2;
		my $quality = $quality1.$quality2;
		
		my $out = '';
		if(exists$hDirectBase{$counter}){
			$out = $hDirectBase{$counter};
		}
		elsif(substr($hDecisionBase{$barcode},-2,2) ne 'UN'){
			$out = $hDecisionBase{$barcode};
		}
		else{
			my $c = '';
			my $T = 0;
			$T += (ord($c)-$OPT_qType) while '' ne ($c=chop$quality);
			if ($T/$barcodeLength >= $MLminQuality){
				$summary{$hDecisionBase{$barcode}}{'untrusted'} --;
				$summaryTotal{'untrusted'} --;
				my %P = ();
				foreach my $idx (@aIndexDual){
					$P{$idx} = &calcML($quality,$barcode,$idx);
				}
				my $firstIdx = (sort{$P{$b}<=>$P{$a}}keys%P)[0];
				if(max(values%P)/sum(values%P) > (1-$MLZ1) ){
					$out = $hCodeToClass{$firstIdx};
					$summary{$hCodeToClass{$firstIdx}}{'trusted'} ++;
					$summaryTotal{'trusted'} ++;
				}else{
					$out = $hCodeToClass{$firstIdx}.'UN';
					$summary{$hCodeToClass{$firstIdx}}{'untrusted'} ++;
					$summaryTotal{'untrusted'} ++;
				}
			}
			else{
				$out = $hDecisionBase{$barcode};
			}
		}

		my $outRead1 = $outhandle{$out."R1"};
		my $outRead2 = $outhandle{$out."R2"};
		my $outIdx1 = $outIdxHandle{$out."I1"};
		my $outIdx2 = $outIdxHandle{$out."I2"};
		print $outRead1 @Fr1;
		print $outRead2 @Fr2;
		print $outIdx1 @Ir1;
		print $outIdx2 @Ir2;
		
		last if eof;
	}
}

close $_ for @inHandle,@idxHandle;
close $outhandle{$_} for keys%outhandle;

print &checkTime()."Output finished.\n";

}

#########################################

#######################
#### FUNCTIONS end ####
#######################

#########################################

__END__
