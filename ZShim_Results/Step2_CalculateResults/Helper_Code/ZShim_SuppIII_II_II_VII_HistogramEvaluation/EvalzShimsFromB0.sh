
do_all() {
	INPUTPATH=/data/p_02098/data/ZSHIM_STUDY/$1                       #raw dicom data
	OUTPUTPATH=/data/pt_02098/ZSHIM_STUDY/Histogram/$1                #output path
	inputDirDicom=`find $INPUTPATH -name 'DICOM'`

	MagnB0Nr=$2
	B0Nr=$3	
	ZrefNr=$4
	
	if [ "$useFM_GRE" == "true" ] ; then
		addStr='_GRE'
		MagnB0Nr=`{ echo "$4 + 1" ; } | bc`
		B0Nr=`{ echo "$4 + 2" ; } | bc`		
	else
		addStr='_CBSFL'
	fi

	ResultsPath=$OUTPUTPATH/Results$addStr'_easy_noGradLimit'

	convert_DICOM_12.sh $inputDirDicom $OUTPUTPATH/B0$addStr $B0Nr

	fslswapdim $OUTPUTPATH/B0$addStr/$B0Nr/timeseries.nii RL PA IS  $OUTPUTPATH/B0$addStr/$B0Nr/timeseries_swapped.nii

	B0map=`find $OUTPUTPATH/B0$addStr -name '*swapped.nii'`
	T2map=$OUTPUTPATH/T2_mask.nii.gz

	if ! [ -d $ResultsPath ] ; then mkdir $ResultsPath ; fi

	{
		echo '.r $IDL_LOCALDIR/PROJECTS/ZSHIMS/EvalzShimsFromB0.pro'
		echo "$ResultsPath"
		echo "$B0map"
		echo "$T2map"
		echo "$useFM_GRE"
	} | idl 
}

export FSLOUTPUTTYPE=NIFTI

useFM_GRE=true
#convert the DICOMs acc to scan numbers

do_all ZS_001 7 9 4
do_all ZS_002 7 9 4
do_all ZS_003 7 9 4
do_all ZS_004 7 9 4
do_all ZS_005 7 9 4
do_all ZS_006 7 9 4
do_all ZS_007 7 9 4
do_all ZS_008 7 9 4
do_all ZS_009 8 10 5
do_all ZS_010 7 9 4
do_all ZS_011 7 9 4
do_all ZS_012 7 9 4
do_all ZS_013 8 10 5
do_all ZS_014 7 9 4
do_all ZS_015 7 9 4
do_all ZS_016 7 9 4
do_all ZS_017 8 10 5
do_all ZS_018 7 9 4
do_all ZS_019 8 10 5
do_all ZS_020 7 9 4
do_all ZS_021 7 9 4
do_all ZS_022 7 9 4
do_all ZS_023 7 9 4
do_all ZS_024 8 10 5
do_all ZS_025 7 9 4
do_all ZS_026 7 9 4
do_all ZS_027 7 9 4
do_all ZS_028 8 10 5
do_all ZS_029 7 9 4
do_all ZS_030 8 10 5
do_all ZS_031 7 9 4
do_all ZS_032 7 9 4
do_all ZS_033 7 9 4
do_all ZS_034 7 9 4
do_all ZS_035 7 9 4
do_all ZS_036 8 10 5
do_all ZS_037 8 10 5
do_all ZS_038 7 9 4
do_all ZS_039 7 9 4
do_all ZS_040 7 9 4
do_all ZS_041 7 9 4
do_all ZS_042 8 10 5
do_all ZS_043 7 9 4
do_all ZS_044 7 9 4
do_all ZS_045 7 9 4
do_all ZS_046 7 9 4
do_all ZS_047 7 9 4
do_all ZS_048 7 9 4
#
