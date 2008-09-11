sub assert_file_same_within($file1,file2){

%status=("status" => $status,
	"line" => $linenr,
	"column" => $collumn,
	"maxerror"=> $maxerror,
	"averageerror"=> $averageerr
	);
return \%status;
}
