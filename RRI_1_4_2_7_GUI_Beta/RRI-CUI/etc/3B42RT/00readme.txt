3B42RT Data Processing

[preparation]
Install "clink" to use bash script (*.sh) on windows command prompt
http://code.google.com/p/clink/

FTP Site for downloading 3B42RT Data
ftp://trmmopen.gsfc.nasa.gov/pub/merged

1) copy header information of your target river basin in ./calc_area/calc_area.txt

2) execute "calc_area.exe" to obtain "out_by_calc_area.txt" file

3) download "3B42RT.20*****.7R2.bin.gz" files from the above FTP site and save them under ./read/infile/

4) "bash unzip.sh" to unzip the downloaded files

5) edit "read_rt_file.sh" file to set extraction range in L4 to L7 (jleft, ibottom, jright, itop) suggested by "out_by_calc_area.txt"

6) "bash read_rt_file.sh" to extract data
note: the extract does not run if output the same output file already exist

7) edit "combine.sh" by setting the extraction range in L4 to L7 (jleft, ibottom, jright, itop) suggested by "out_by_calc_area.txt", and set output file name on L9. Also edit L14, L17 and L20 to indicate which year, month and day of the data should be processed.

8) "bash ./combine.sh" to combine all rainfall files to create the RRI input
