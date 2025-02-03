my $name = $ARGV[1];
        open FILE, "+<$name"
          or die "Cannot open $name: $!";
        flock FILE, 2;
 	local $/=undef;
        $data =  <FILE>;
	$data=~s/(^\s*subroutine.+?$)/$1\n$ARGV[0]/gmi;
        seek FILE, 0, 0;
        truncate FILE, 0;
        print FILE $data;
 	print $data;
        close FILE;
