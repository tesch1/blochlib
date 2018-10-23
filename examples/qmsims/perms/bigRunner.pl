#!/usr/bin/perl 

#<?php

#----------------------------------
#several data directories
$dataDir="data/";  #place to store output data
$tensorsDir="traceMats/"; #the binary tensor files
$permCyDir="permutationTrains/"; #the permutation file dir
$confDir="config_files/"; #place to store config files

$tmpDir="tmp/"; # tmp dir 

#----------------------------------
#base file names
$pCybase="pCy_";  #the permutation file base
$nameBase="n"; #base name for the names ouput
$seqBase="s"; #base names for the sequence output

#----------------------------------
#----------------------------------
#the spin parameter section variable

$c13Pars="
	T 13C 0
	T 13C 1 
	D 2146 0 1 0 0 0 
	C 1254 12345 0 0 
	C -1544 8552 0 1\n";

%spinList=(
	0=>"numspin 2".$c13Pars,
	
	
	1=>"numspin 3".$c13Pars."
	T 1H 2
	D 4506 0 2 0 0 0 
	D 7564 1 2 0 0 0 ",
	
	2=>"numspin 4".$c13Pars."
	T 1H 2
	T 1H 3
	D 4506 0 2 0 0 0 
	D 7564 1 2 0 0 0 
	D 15546 2 3 0 0 0 
	D 2150 0 3 0 0 0 
	D 4562 1 3 0 0 0 "
);
	
$num13C="2";
$num1H="0";


#----------------------------------
#parameter variables

#rotor speed and angle
$wr=5000;
$rotor=54.7356;


#the pregenerated binary tensor file
# generated from the 'tensorgen' program
#can be more then one file...just separate file names with a comma
#available tensor files
@tensorFiles=("cart2.bin", "cart3.bin", "sphere2.bin", "sphere3.bin", "hamil.bin");
$tensorfile=$tensorsDir.$tensorFiles[0].",".$tensorsDir.$tensorFiles[2];

# largest time for a propogator step 
$maxtimestep=1e-6;

#powder average type	
#powder ../../../cyrstals/rep256.cry
$powder="ZCW_3_1154";

#for non-iquid powder average types
$thetasteps =616;	
$phisteps=377;


# if this is '-1' all will be calculated
# if this is '-2' NONE will be calculated
# else if will calculate the trace for that train only
$traceForIndex=-1;

#matlab file where data is saved
# naming is as  "<recouple type>_p<permution length>_n<napps>_<rotorspeed>_c<num 13C>_h<num1H>"
# so a sinlge post-C7 sequnce applied 10 time would be
# spinning at 10 Khz
# "Cp71_p1_n10_10k.mat"
#matout Cp71_p1_n32_5k_csa2.mat
#matout CpCppBar71_p1_n8_5k_csa2_4c_a.mat
#matout Cp71_p1_n20_5k_csa2.mat
$matout=""; #to be set by the runner

#output file name for the seqeunce names
# should be 'blaaa'.m a simple matlab script file
$seqout=""; #to be set by the runner

#this is the name of the matlab file
# that makes an array of the names of the tenor types
# so that you can load it as a function
$nameout=""; #to be set by the runner

#if you deisre to collect an fid of a certain train
# set to '-1' if no fid desired
$fidForIndex=-1;


#debug level (info dumping) level for the program
$debuglevel=0;
$progress=0;

#inary output file and options for generation
# can be: ReadAndRun, GenerateAndSave, GenerateAndRun
$generationOps="ReadAndRun";


$trainFile=""; #to be set by the runner

#----------------------------------
#pulses variables

$pulsespin="13C";
$baseSym=7;
$factorSym=1;

#number of rotor periods the sequence should be applied over
$wrfactor=2;

#number of applications of the sequence
$napps=1;

#subunit ids define what chars will be used as units below

@subUnits=('o', 'O', 'w', 'W');
$subUnitIds=""; #to be set by runner

#types associate a 'rtype' with a subunit id
# valid types are: R, Rp, RBar, RpBar, 
# C, CpBar, CBar, Cpp, CppBar
%Types=(
	'o'=>"Cp",
	'O'=>"CpBar",
	'w'=>"Cpp",
	'W'=>"CppBar"
);

$typeIds=""; #to be set by runner

%UnitsList=(
	
	1=>{ "a" => ["o",1] },
	
	2 =>{
		"a" => ["o, O", 1],
		"d" => ["w, W",1],
		"g" => ["oO, wW",2],		
		"h" => ["wO, oW",2],		
		"i" => ["ow ,OW",2],
		"j" => ["oOOo, wWWw",4],
		"k" => ["oOWw, OowW",4],
		"l" => ["oOOo, WwwW",4]
		
		},	
	4 =>{
		"a" => [ "o, O, w, W",1],
		"b" =>  ["oOOo, wWWw, Ww, Oo",6]
		}
);

$units=""; #to be set by runner

#the 'longest' list to permute the
$permutations=1; #to be set by runner

%masBy=(
	1=>[2, 4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48],
	2=>[2,4,8,12,16,20], #those with 2xN trainfiles
	4=>[4,8,12]	#those with 4xN trainfiles
);




#----------------------------------

#the out put file opener
sub newopen {
	my $path = shift;
	local  *FH;  # not my!
	open   (FH, $path) 	    or  return undef;
	return *FH;
}

#master output file
$FOUT = *STDOUT;

#makes the two strings 
# subUnitIds o...
# typeIds Cp...
# based on what chars are in the 'units'
sub genSubStrs{
	my $pL = shift;
	local $outStr="", $item ;
	local $subUss="subUnitIds ";
	local $typeUss="typeIds ";
	local @coo, $item;
	$_=$pL;
	foreach $item (keys %Types){
		if(m/$item/ ){
			$subUss.=$item.",";
			$typeUss.= $Types{$item} .",";
		}		
	}
	chop($subUss);
	chop($typeUss);
	$outStr=$subUss."\n\t".$typeUss;	
	return $outStr;
}
	


#the master runner

#first we need to get all the 'normal' Cp experiments done
# essentially 
#here 
#permutations=1,  subUnitIds='o', napps {2 or 4}x{by2x or by4x},
#typeIds=Cp

#first dump the spininfo
$BaseMatOuts=[];
$BaseNames=[];
$BaseSeqs=[];
$NappsMat=[];

$maxSpins=1; $onSpin=0;
$maxOrder=1; $onOrder=0;

$i=-1;

#open the PBS que file

#this tells me when to start and end 'printing'
# the jobs (useful for when you have to kill
# a large config file, but many have already run...)
$startPBSct=0;
$endPBSct="all";
$curPBSct=0;


	
$pbsCT=0;

foreach $num1H (keys %spinList){  #spin system tests
#	$onSpin=$onSpin+1;
#	if($onSpin>$maxSpins){ exit;	}
	
	#Save the matlab files of names for p=2 and 4
	local *matNames,*matNamesNorms;
	open(matNames, '>',$dataDir.'filenames_h'.$num1H.'.m');

	#Save the matlab files of names for p=1
	open(matNamesNorms, '>',$dataDir.'filenames_norm_h'.$num1H.'.m');
	
$i=-1;

foreach $basV ( keys %masBy ){  #the 'order number' (2 or 4)
	
	close(PBS);
	$pbsfile="pbsfiles/pbs_".$basV."_".$flags."_".$num1H.".sh";
	open(PBS, '>', "".$pbsfile);
	print {PBS} "
	#PBS -o waugh.cchem.berkeley.edu:/usr/people/magneto/code/SIMS/perms/tmp/
	#PBS -l ncpus=4
	#PBS -A magneto
	#PBS -q longRun\@waugh.cchem.berkeley.edu
	#PBS -N MIghtyRec
	#PBS -j oe
	
	cd /usr/people/magneto/code/SIMS/perms
\n";
#	$onOrder=$onOrder+1;
#	if($onOrder>$maxOrder){ exit;	}
	

foreach $item ( 0..$#{$masBy{$basV}  }){ #the number of permutations
	$item= $masBy{$basV}[$item];

for $flags ( keys %{ $UnitsList{$basV} } ){ #the Units list
	
	$sunit=$UnitsList{$basV}{$flags}[0];
	$mulunit=$UnitsList{$basV}{$flags}[1];
	$typesAndSubU=genSubStrs($sunit);
	
	$i=$i+1;	
	
	if($basV==1){
		$napps=($item);
		$permutations=1;
	}else{
		$napps=1;
		$permutations=$item;
	}
	$tail="Cp".$baseSym.$factorSym."_p".$basV."x".$item.
	      "_n".$napps."_w".$wr."_c".$num13C.
	      "_h".$num1H."_f".$flags;
	
	$matout=$dataDir.$tail.".mat";
	
	$nameout=$dataDir.$nameBase."_".$tail.".m";
	$seqout=$dataDir.$seqBase."_".$tail.".m";
	
	$BaseMatOuts->[$i]=$matout;
	$BaseNames->[$i]=$nameBase."_".$tail;
	$BaseSeqs->[$i]=$seqBase."_".$tail;
	if($basV==1){
		$NappsMat->[$i]=$napps;
	}else{
		$NappsMat->[$i]=($item)*$mulunit;
	}
	
	
	$trainFile=$permCyDir.$pCybase.$basV."x".$item.".bn";
	
	$confFile=$confDir."rc_".$basV."x".$item."_".$flags."_h".$num1H.".sim";
	open(fout, '>', $confFile);
	
	print {fout} "
spins{ 
	".$spinList{$num1H}."
}

params{
#rotor
	wr $wr
	rotor $rotor
	
#the pregenerated binary tensor file
# generated from the 'tensorgen' program
#can be more then one file...just separate file names with a comma
	tensorfile $tensorfile

#powder
	powder $powder
	thetasteps $thetasteps
	phisteps $phisteps

#integration max time
	maxtimestep $maxtimestep

#index to trace
	traceForIndex	$traceForIndex

#fid for index
	fidForIndex $fidForIndex

#output levels
	debuglevel $debuglevel
	progress $progress

#data output names
	matout $matout
	nameout $nameout
	seqout $seqout

#trains file
	trainFile $trainFile

#what to do
	generationOps $generationOps
}

pulses{
	pulsespin $pulsespin

	baseSym $baseSym
	factorSym $factorSym

#number of rotor periods the sequence should be applied over
	wrfactor $wrfactor

	napps $napps

	$typesAndSubU

	units $sunit

	permutations $permutations
}

";
	close(fout);

	$curPBSct=$curPBSct+1;
	print {STDERR} "$curPBSct | $confFile\n";
	if($curPBSct>=$startPBSct and 
	  ($curPBSct<$endPBSct or $endPBSct=="all"))
	{
		$pbsCT=$pbsCT+1;
		#print {PBS} "echo $confFile\n";
		if($pbsCT<4){
			print {PBS} "rectrain $confFile &\n";		
		}else{
			print {PBS} "rectrain $confFile \n";
			$pbsCT=0;
		}
	}
		
	

} #end arr array


} #end the Units list

	#---- BEGIN MAT FILES NAMES OUT----
	# dump out the files names we will record
	$curFile=matNames;
	
	$totpath=`pwd`; chop($totpath);

	if($basV!=4){ $curFile=NULL;	}
	if($basV==1){ $curFile=matNamesNorms;	}

	print {$curFile} "path(path, '$totpath/$dataDir');\n";
	print {$curFile} "napps=[";
	for($j=0;$j<=$i;$j=$j+1){
		if($j==($i)){	
			print {$curFile} "".$NappsMat->[$j]."];\n\n";
		}else{
			print {$curFile} "".$NappsMat->[$j].",\n";
		}
	}

	print {$curFile} "namefiles=cell($i+1,1); \n";
	print {$curFile} "namefiles={";
	for($j=0;$j<=$i;$j=$j+1){
		if($j==($i)){	
			print {$curFile} "'".$BaseNames->[$j]."'};\n\n";
		}else{
			print {$curFile} "'".$BaseNames->[$j]."',\n";
		}
	}

	print {$curFile} "seqfiles=cell($i+1,1); \n";
	print {$curFile} "seqfiles={";
	for($j=0;$j<=$i;$j=$j+1){
		if($j==($i)){	
			print {$curFile} "'".$BaseSeqs->[$j]."'};\n\n";
		}else{
			print {$curFile} "'".$BaseSeqs->[$j]."',\n";
		}
	}

	print {$curFile} "matfiles=cell($i+1,1); \n";
	print {$curFile} "matfiles={";
	for($j=0;$j<=$i;$j=$j+1){
		if($j==($i)){	
			print {$curFile} "'".$BaseMatOuts->[$j]."'};\n\n";
		}else{
			print {$curFile} "'".$BaseMatOuts->[$j]."',\n";
		}
	}

	if($basV==1){ 
		$BaseMatOuts=[];
		$BaseNames=[];
		$BaseSeqs=[];
		$NappsMat=[];
		$i=-1;
	}

} #end list hash

#---- END MAT FILES NAMES OUT----

} #end spinsSys hash

#spawn the job
#system("qsub $pbsfile");



#?>
