sub assert_file_same_within
{
	$file1=$_[0],file2=$_[1]
	open FILE1 $file1
	open FILE2 $file2
	while(<FILE1>)
	{
		$linefile1=$_
		$linefile2=<fILE2>
		while($linefile1=~m/$number/)
		{
			push(@numbers1,$_)	
			$linefile2=~m/$number/
			push(@numbers2,$_)
		}
	}

if(max(@numbers1-@numbers2)>$tol)
{
 $status=0
 $maxerror=max(@numbers1-@numbers2)
 }
%test=("status" => $status,
	"line" => $linenr,
	"column" => $collumn, 
	"maxerror"=> $maxerror,
	"averageerror"=> $averageerr
	);
return \%status;
}
