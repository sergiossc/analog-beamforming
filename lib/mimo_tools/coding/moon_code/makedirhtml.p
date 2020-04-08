#!/usr/bin/perl

# create an html file with the diretocty listing
# makedirhtml directory 
$dir = shift;
$outputfile = shift;

$outputfile = "index.html";
unlink($outputfile);

open(OUTF,">".$outputfile) || die "Can't open $outputfile\n";

# Print header stuff
print OUTF "<!doctype html public>\n";
print OUTF "<meta name=\"author\" content=\"Todd K. Moon\">\n";
print OUTF "<meta name=\"copyright\" content=\"2003 Todd K. Moon\">\n";
print OUTF "<html>\n";
print OUTF "<head>\n";
print OUTF "<title> Index of files </title>\n";
print OUTF "</head>\n";
print OUTF "<body>\n";
print OUTF "<h1> Index of files </h1>\n";

print OUTF "<ul>\n";


foreach $dl (`ls *`) {
	chop $dl;
	print $dl,"\n";
    
   if(!($dl =~ m-^/-)) {       # it is a not directory
	   $fname = $dl;
   }
	print $fname,"\n";
	print OUTF "<li> <a href=\"$fname\"> $fname </a> </li>\n";
}

print OUTF "</ul>\n";
print OUTF "</body>\n";
print OUTF "</html>\n";

	   
close OUTF;
