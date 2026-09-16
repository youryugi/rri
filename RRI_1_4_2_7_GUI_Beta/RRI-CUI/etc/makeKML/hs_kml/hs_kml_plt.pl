use strict;
use warnings;


foreach my $id (1..96){

 my $num1 = sprintf("%06d",$id);
 my $num2 = sprintf("%03d",$id);

 my $argment = "-e \"s/number/$num1/g\" ";
 system("sed $argment hs_kml.plt > hs_temp.plt");

 $argment = "-e \"s/tstep/$num2/g\" ";
 system("sed $argment hs_temp.plt > hs_$num2.plt");

 print "hs_$num1.plt\n";

 system("gnuplot hs_$num2.plt");

 #system("rm hs_$num1.plt");

}



#for /l %%i in (1, 1, 152) do (
# set num=00000%%i
# set num2=!num:~-6!
# sed -e "s/number/!num2!/g" hs.plt > hs_temp.plt

# set num3=!num:~-3!
# sed -e "s/tstep/!num3!/g" hs_temp.plt > hs_!num2!.plt

# echo hs_!num2!.plt
# wgnuplot hs_!num2!.plt
# del hs_!num2!.plt
# del hs_temp.plt

#)
