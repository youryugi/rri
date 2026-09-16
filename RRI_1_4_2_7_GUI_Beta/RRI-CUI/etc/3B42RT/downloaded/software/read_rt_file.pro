;---------------------------------------------------------------
;
;	filename:	read_rt_file.pro
;
;	purpose:
;
;	This IDL source file contains the read_3B4XRT
;	procedure, which reads all the data in a 3B40RT,
;	3B41RT, or 3B42RT binary data file.  Recall that the data
;	files vary in latitudinal extent depending on the file
;	type specified.
;
;	The read_3B4XRT procedure swaps the bytes of the data arrays 
;	if necessary.
;
;	This source file also contains the read_3B4XRT_header
;	procedure which reads just the header of one of these
;	binary files.
;
;	Instructions:
;
;	Type "idl" at the UNIX command prompt to start IDL.
;	Type ".run read_rt_file.pro" to compile the procedures.
;	Type "read_3B4XRT, file, data" where "file" is the
;	name of the file to be read and "data" is the data
;	structure that will contain the data read from the file
;	in two fields, "data.header" and "data.precip".
;
;	change log:
;
;	O. Kelley		02/20/2002
;	G.J. Huffman		06/20/2003 Correct header size
;       G.J. Huffman            05/14/2007 Add swapping for i386
;                                          machines (S. Nesbitt)
;	G.J. Huffman		02/25/2008 Documentation
;	G.J. Huffman		01/30/2009 Update for Feb. 2009
;					   TMPA-RT release
;---------------------------------------------------------------


	pro read_3B4XRT_alg_ID, file, alg_ID
;-----------------------------------------
;	The "alg_ID" metadata is algorithm_ID=xxxxxx, where xxxxxx 
;	is the algorithm identifier.  Nence the "+13" in strmid.
;-----------------------------------------
	num_lon		= 1440
	header		= bytarr( num_lon*2 )

	close,1
	openr,1, file
	readu,1, header
	close,1

	header		= strtrim( string(header), 2 )
	alg_ID		= strpos( header, 'algorithm_ID=' )
	alg_ID		= strmid( header, alg_ID+13, 6 )

	end

	pro create_3B4XRT_struct, alg_ID, data
;-----------------------------------------
;	Set up the data structure for output.
;-----------------------------------------
	num_lon		= 1440
	num_lat_90	=  720
	num_lat_60	=  480

	case alg_ID of
	  '3B40RT': data = { header:	bytarr(num_lon*2), $
	    precip:		intarr(num_lon,num_lat_90), $
	    precip_error:	intarr(num_lon,num_lat_90), $
	    total_pixels:	bytarr(num_lon,num_lat_90), $
	    ambiguous_pixels:	bytarr(num_lon,num_lat_90), $
	    rain_pixels:	bytarr(num_lon,num_lat_90), $
	    source_of_estimate:	bytarr(num_lon,num_lat_90)  $
	  }

	  '3B41RT': data = { header:	bytarr(num_lon*2), $
	    precip:		intarr(num_lon,num_lat_60), $
	    precip_error:	intarr(num_lon,num_lat_60), $
	    total_pixels:	bytarr(num_lon,num_lat_60)  $
	  }

	  '3B42RT': data = { header:	bytarr(num_lon*2), $
	    precip:		intarr(num_lon,num_lat_60), $
	    precip_error:	intarr(num_lon,num_lat_60), $
	    source_of_estimate:	bytarr(num_lon,num_lat_60), $
	    precip_uncal:	intarr(num_lon,num_lat_60)  $
	  }

	  else: print, 'create_3B4XRT_struct: unknown alg_ID... ', alg_ID
	endcase

	end

	pro byte_swap_3B4XRT, data
;-----------------------------------------
;	Swap bytes if needed.
;-----------------------------------------
;
; ----- Determine byte order of system and of file (keyed to "byte_order=
;	big_endian").
;
         arch            = strlowcase( !version.arch )

         if (arch eq 'x86') OR (arch eq 'alpha') OR (arch eq 'i386') $
           then system_byte_order       = 'little_endian' $
           else system_byte_order       = 'big_endian'

	if strpos( string(data.header), 'big_endian') ne -1 $
	  then file_byte_order		= 'big_endian' $
	  else file_byte_order		= 'little_endian'
;
; ----- if necessary, swap the bytes of the 2-byte variables
;
	if system_byte_order ne file_byte_order then begin
	  print, 'byte_swap_3B4XRT: warning: swapping bytes...'
	  data.precip		= swap_endian( data.precip )
	  data.precip_error	= swap_endian( data.precip_error )
	endif

	end

	pro read_3B4XRT, file, data
;-----------------------------------------
;	The main procedure; get the algorithm identifier, create the 
;	data structure, read both header and data, and swap bytes if 
;	needed.
;-----------------------------------------

	read_3B4XRT_alg_ID, file, alg_ID
	create_3B4XRT_struct, alg_ID, data

	close,1
	openr,1, file
	readu,1, data
	close,1

	byte_swap_3B4XRT, data

	end

	pro read_3B4XRT_header, file, header
;-----------------------------------------
;	The procedure to just read the header.
;-----------------------------------------
	num_lon		= 1440
	header		= bytarr( num_lon*2 )

	close,1
	openr,1, file
	readu,1, header
	close,1

	header		= str_sep( strtrim(string(header),2), ' ' )

	end

;-----------------------------------------
;	Help.
;-----------------------------------------

	print, ' '
	print, '  usage: '
	print, '    read_3B4XRT, file, data'
	print, '    read_3B4XRT_header, file, header'
	print, ' '

	end
