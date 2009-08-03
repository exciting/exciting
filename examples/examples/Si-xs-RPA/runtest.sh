#! /bin/sh

################################################################################

# files to compare for the test

files="\
	*.in \
	EFERMI.OUT EIGVAL.OUT EQATOMS.OUT EVALCORE.OUT FERMIDOS.OUT GEOMETRY.OUT \
	IADIST.OUT INFO.OUT KPOINTS.OUT LATTICE.OUT LINENGY.OUT RMSDVEFF.OUT \
	SYMCRYS.OUT SYMLAT.OUT SYMSITE.OUT TOTENERGY.OUT \
	\
	BAND*.OUT \
	\
	PARAMS*.OUT INFOXS.OUT TIME.OUT \
	SYMGENR.OUT SYMINV.OUT SYMMULT.OUT SYMMULT_TABLE.OUT SYMT2.OUT \
	GQPOINTS_*.OUT KMAPKQ_*.OUT QPOINTS*.OUT \
	EIGVAL_*.OUT EFERMI_*.OUT \
	FXC_BSE_HEAD*.OUT \
	BSEDIAG.OUT DIELTENS*_*.OUT \
	EXCITON_*.OUT EPSILON_*.OUT LOSS_*.OUT SIGMA_*.OUT SUMRULES_*.OUT \
"
# The following files are larger than the others in the list above.
# They can be commented out
files_large="\
	SCREEN_*.OUT W_SCREEN_*.OUT EXCLI_ASC.OUT SCCLI_ASC.OUT
"
files="$files $files_large"

################################################################################


rundir="_current"
diffdir="_diff.reference_current"

echo
echo "Running test..."
echo

\cp -f exciting.in_ exciting.in

../../bin/excitingser

\rm -fr $rundir
\rm -fr $diffdir
\mkdir $rundir
\mkdir $diffdir

for file in $files; do
	f=$file"_"
	\cp -f $file $rundir/$f
	diff -B -w _reference/$f $rundir/$f > $diffdir/$file.diff
done

echo "done."
echo

