# Adapted from libtiff v4.0.9

#   Copyright (c) 1988-1997 Sam Leffler
#   Copyright (c) 1991-1997 Silicon Graphics, Inc.

#   Permission to use, copy, modify, distribute, and sell this software and
#   its documentation for any purpose is hereby granted without fee, provided
#   that (i) the above copyright notices and this permission notice appear in
#   all copies of the software and related documentation, and (ii) the names of
#   Sam Leffler and Silicon Graphics may not be used in any advertising or
#   publicity relating to the software without the specific, prior written
#   permission of Sam Leffler and Silicon Graphics.

#   THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
#   WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

#   IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
#   ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
#   OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
#   WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF
#   LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
#   OF THIS SOFTWARE.


#   NB: In the comments below,
#    - items marked with a + are obsoleted by revision 5.0,
#    - items marked with a ! are introduced in revision 6.0.
#    - items marked with a % are introduced post revision 6.0.
#    - items marked with a $ are obsoleted by revision 6.0.
#    - items marked with a & are introduced by Adobe DNG specification.



"""
	$(TYPEDEF)

List of many common named TIFF Tags. This is not an exhaustive list but should
cover most cases.

$(FIELDS)

"""
@enum TiffTag begin
	SUBFILETYPE = 254 			 #  subfile data descriptor
#define	    FILETYPE_REDUCEDIMAGE	0x1	/* reduced resolution version */
#define	    FILETYPE_PAGE		0x2	/* one page of many */
#define	    FILETYPE_MASK		0x4	/* transparency mask */
	OSUBFILETYPE = 255 			 #  +kind of data in subfile
#define	    OFILETYPE_IMAGE		1	/* full resolution image data */
#define	    OFILETYPE_REDUCEDIMAGE	2	/* reduced size image data */
#define	    OFILETYPE_PAGE		3	/* one page of many */
	IMAGEWIDTH = 256 			 #  image width in pixels
	IMAGELENGTH = 257 			 #  image height in pixels
	BITSPERSAMPLE = 258 			 #  bits per channel (sample)
	COMPRESSION = 259 			 #  data compression technique
	PHOTOMETRIC = 262 			 #  photometric interpretation
	THRESHHOLDING = 263 			 #  +thresholding used on data
#define	    THRESHHOLD_BILEVEL		1	/* b&w art scan */
#define	    THRESHHOLD_HALFTONE		2	/* or dithered scan */
#define	    THRESHHOLD_ERRORDIFFUSE	3	/* usually floyd-steinberg */
	CELLWIDTH = 264 			 #  +dithering matrix width
	CELLLENGTH = 265 			 #  +dithering matrix height
	FILLORDER = 266 			 #  data order within a byte
#define	    FILLORDER_MSB2LSB		1	/* most significant -> least */
#define	    FILLORDER_LSB2MSB		2	/* least significant -> most */
	DOCUMENTNAME = 269 			 #  name of doc. image is from
	IMAGEDESCRIPTION = 270 			 #  info about image
	MAKE = 271 			 #  scanner manufacturer name
	MODEL = 272 			 #  scanner model name/number
	STRIPOFFSETS = 273 			 #  offsets to data strips
	ORIENTATION = 274 			 #  +image orientation
#define	    ORIENTATION_TOPLEFT		1	/* row 0 top, col 0 lhs */
#define	    ORIENTATION_TOPRIGHT	2	/* row 0 top, col 0 rhs */
#define	    ORIENTATION_BOTRIGHT	3	/* row 0 bottom, col 0 rhs */
#define	    ORIENTATION_BOTLEFT		4	/* row 0 bottom, col 0 lhs */
#define	    ORIENTATION_LEFTTOP		5	/* row 0 lhs, col 0 top */
#define	    ORIENTATION_RIGHTTOP	6	/* row 0 rhs, col 0 top */
#define	    ORIENTATION_RIGHTBOT	7	/* row 0 rhs, col 0 bottom */
#define	    ORIENTATION_LEFTBOT		8	/* row 0 lhs, col 0 bottom */
	SAMPLESPERPIXEL = 277 			 #  samples per pixel
	ROWSPERSTRIP = 278 			 #  rows per strip of data
	STRIPBYTECOUNTS = 279 			 #  bytes counts for strips
	MINSAMPLEVALUE = 280 			 #  +minimum sample value
	MAXSAMPLEVALUE = 281 			 #  +maximum sample value
	XRESOLUTION = 282 			 #  pixels/resolution in x
	YRESOLUTION = 283 			 #  pixels/resolution in y
	PLANARCONFIG = 284 			 #  storage organization
#define	    PLANARCONFIG_CONTIG		1	/* single image plane */
#define	    PLANARCONFIG_SEPARATE	2	/* separate planes of data */
	PAGENAME = 285 			 #  page name image is from
	XPOSITION = 286 			 #  x page offset of image lhs
	YPOSITION = 287 			 #  y page offset of image lhs
	FREEOFFSETS = 288 			 #  +byte offset to free block
	FREEBYTECOUNTS = 289 			 #  +sizes of free blocks
	GRAYRESPONSEUNIT = 290 			 #  $gray scale curve accuracy
#define	    GRAYRESPONSEUNIT_10S	1	/* tenths of a unit */
#define	    GRAYRESPONSEUNIT_100S	2	/* hundredths of a unit */
#define	    GRAYRESPONSEUNIT_1000S	3	/* thousandths of a unit */
#define	    GRAYRESPONSEUNIT_10000S	4	/* ten-thousandths of a unit */
#define	    GRAYRESPONSEUNIT_100000S	5	/* hundred-thousandths */
	GRAYRESPONSECURVE = 291 			 #  $gray scale response curve
	T4OPTIONS = 292 			 #  TIFF 6.0 proper name alias
#define	    GROUP3OPT_2DENCODING	0x1	/* 2-dimensional coding */
#define	    GROUP3OPT_UNCOMPRESSED	0x2	/* data not compressed */
#define	    GROUP3OPT_FILLBITS		0x4	/* fill to byte boundary */
	T6OPTIONS = 293 			 #  TIFF 6.0 proper name
#define	    GROUP4OPT_UNCOMPRESSED	0x2	/* data not compressed */
	RESOLUTIONUNIT = 296 			 #  units of resolutions
#define	    RESUNIT_NONE		1	/* no meaningful units */
#define	    RESUNIT_INCH		2	/* english */
#define	    RESUNIT_CENTIMETER		3	/* metric */
	PAGENUMBER = 297 			 #  page numbers of multi-page
	COLORRESPONSEUNIT = 300 			 #  $color curve accuracy
#define	    COLORRESPONSEUNIT_10S	1	/* tenths of a unit */
#define	    COLORRESPONSEUNIT_100S	2	/* hundredths of a unit */
#define	    COLORRESPONSEUNIT_1000S	3	/* thousandths of a unit */
#define	    COLORRESPONSEUNIT_10000S	4	/* ten-thousandths of a unit */
#define	    COLORRESPONSEUNIT_100000S	5	/* hundred-thousandths */
	TRANSFERFUNCTION = 301 			 #  !colorimetry info
	SOFTWARE = 305 			 #  name & release
	DATETIME = 306 			 #  creation date and time
	ARTIST = 315 			 #  creator of image
	HOSTCOMPUTER = 316 			 #  machine where created
	PREDICTOR = 317 			 #  prediction scheme w/ LZW
#define     PREDICTOR_NONE		1	/* no prediction scheme used */
#define     PREDICTOR_HORIZONTAL	2	/* horizontal differencing */
#define     PREDICTOR_FLOATINGPOINT	3	/* floating point predictor */
	WHITEPOINT = 318 			 #  image white point
	PRIMARYCHROMATICITIES = 319 			 #  !primary chromaticities
	COLORMAP = 320 			 #  RGB map for palette image
	HALFTONEHINTS = 321 			 #  !highlight+shadow info
	TILEWIDTH = 322 			 #  !tile width in pixels
	TILELENGTH = 323 			 #  !tile height in pixels
	TILEOFFSETS = 324 			 #  !offsets to data tiles
	TILEBYTECOUNTS = 325 			 #  !byte counts for tiles
	BADFAXLINES = 326 			 #  lines w/ wrong pixel count
	CLEANFAXDATA = 327 			 #  regenerated line info
#define	    CLEANFAXDATA_CLEAN		0	/* no errors detected */
#define	    CLEANFAXDATA_REGENERATED	1	/* receiver regenerated lines */
#define	    CLEANFAXDATA_UNCLEAN	2	/* uncorrected errors exist */
	CONSECUTIVEBADFAXLINES = 328 			 #  max consecutive bad lines
	SUBIFD = 330 			 #  subimage descriptors
	INKSET = 332 			 #  !inks in separated image
#define	    INKSET_CMYK			1	/* !cyan-magenta-yellow-black color */
#define	    INKSET_MULTIINK		2	/* !multi-ink or hi-fi color */
	INKNAMES = 333 			 #  !ascii names of inks
	NUMBEROFINKS = 334 			 #  !number of inks
	DOTRANGE = 336 			 #  !0% and 100% dot codes
	TARGETPRINTER = 337 			 #  !separation target
	EXTRASAMPLES = 338 			 #  !info about extra samples
#define	    EXTRASAMPLE_UNSPECIFIED	0	/* !unspecified data */
#define	    EXTRASAMPLE_ASSOCALPHA	1	/* !associated alpha data */
#define	    EXTRASAMPLE_UNASSALPHA	2	/* !unassociated alpha data */
	SAMPLEFORMAT = 339 			 #  !data sample format
#define	    SAMPLEFORMAT_UINT		1	/* !unsigned integer data */
#define	    SAMPLEFORMAT_INT		2	/* !signed integer data */
#define	    SAMPLEFORMAT_IEEEFP		3	/* !IEEE floating point data */
#define	    SAMPLEFORMAT_VOID		4	/* !untyped data */
#define	    SAMPLEFORMAT_COMPLEXINT	5	/* !complex signed int */
#define	    SAMPLEFORMAT_COMPLEXIEEEFP	6	/* !complex ieee floating */
	SMINSAMPLEVALUE = 340 			 #  !variable MinSampleValue
	SMAXSAMPLEVALUE = 341 			 #  !variable MaxSampleValue
    CLIPPATH = 343           # %ClipPath [Adobe TIFF technote 2]
    XCLIPPATHUNITS = 344     # %XClipPathUnits [Adobe TIFF technote 2]
    YCLIPPATHUNITS = 345	 # %YClipPathUnits [Adobe TIFF technote 2]
    INDEXED	= 346        # %Indexed [Adobe TIFF Technote 3]
	JPEGTABLES = 347 			 #  %JPEG table stream
	OPIPROXY = 351 			 #  %OPI Proxy [Adobe TIFF technote]
# Tags 400-435 are from the TIFF/FX spec
	GLOBALPARAMETERSIFD = 400 			 #  !
	PROFILETYPE = 401 			 #  !
#define     PROFILETYPE_UNSPECIFIED	0	/* ! */
#define     PROFILETYPE_G3_FAX		1	/* ! */
	FAXPROFILE = 402 			 #  !
#define     FAXPROFILE_S			1	/* !TIFF/FX FAX profile S */
#define     FAXPROFILE_F			2	/* !TIFF/FX FAX profile F */
#define     FAXPROFILE_J			3	/* !TIFF/FX FAX profile J */
#define     FAXPROFILE_C			4	/* !TIFF/FX FAX profile C */
#define     FAXPROFILE_L			5	/* !TIFF/FX FAX profile L */
#define     FAXPROFILE_M			6	/* !TIFF/FX FAX profile LM */
	CODINGMETHODS = 403 			 #  !TIFF/FX coding methods
#define     CODINGMETHODS_T4_1D		(1 << 1)	/* !T.4 1D */
#define     CODINGMETHODS_T4_2D		(1 << 2)	/* !T.4 2D */
#define     CODINGMETHODS_T6		(1 << 3)	/* !T.6 */
#define     CODINGMETHODS_T85 		(1 << 4)	/* !T.85 JBIG */
#define     CODINGMETHODS_T42 		(1 << 5)	/* !T.42 JPEG */
#define     CODINGMETHODS_T43		(1 << 6)	/* !T.43 colour by layered JBIG */
	VERSIONYEAR = 404 			 #  !TIFF/FX version year
	MODENUMBER = 405 			 #  !TIFF/FX mode number
	DECODE = 433 			 #  !TIFF/FX decode
	IMAGEBASECOLOR = 434 			 #  !TIFF/FX image base colour
	T82OPTIONS = 435 			 #  !TIFF/FX T.82 options

# Tags 512-521 are obsoleted by Technical Note #2 which specifies a
# revised JPEG-in-TIFF scheme.

	JPEGPROC = 512 			 #  !JPEG processing algorithm
#define	    JPEGPROC_BASELINE		1	/* !baseline sequential */
#define	    JPEGPROC_LOSSLESS		14	/* !Huffman coded lossless */
	JPEGIFOFFSET = 513 			 #  !pointer to SOI marker
	JPEGIFBYTECOUNT = 514 			 #  !JFIF stream length
	JPEGRESTARTINTERVAL = 515 			 #  !restart interval length
	JPEGLOSSLESSPREDICTORS = 517 			 #  !lossless proc predictor
	JPEGPOINTTRANSFORM = 518 			 #  !lossless point transform
	JPEGQTABLES = 519 			 #  !Q matrix offsets
	JPEGDCTABLES = 520 			 #  !DCT table offsets
	JPEGACTABLES = 521 			 #  !AC coefficient offsets
	YCBCRCOEFFICIENTS = 529 			 #  !RGB -> YCbCr transform
	YCBCRSUBSAMPLING = 530 			 #  !YCbCr subsampling factors
	YCBCRPOSITIONING = 531 			 #  !subsample positioning
#define	    YCBCRPOSITION_CENTERED	1	/* !as in PostScript Level 2 */
#define	    YCBCRPOSITION_COSITED	2	/* !as in CCIR 601-1 */
	REFERENCEBLACKWHITE = 532 			 #  !colorimetry info
	STRIPROWCOUNTS = 559 			 #  !TIFF/FX strip row counts
	XMLPACKET = 700	         # %XML packet [Adobe XMP Specification, January 2004]
    OPIIMAGEID = 32781	     # %OPI ImageID [Adobe TIFF technote]
# tags 32952-32956 are private tags registered to Island Graphics
	REFPTS = 32953 			 #  image reference points
	REGIONTACKPOINT = 32954 			 #  region-xform tack point
	REGIONWARPCORNERS = 32955 			 #  warp quadrilateral
	REGIONAFFINE = 32956 			 #  affine transformation mat
# tags 32995-32999 are private tags registered to SGI
	MATTEING = 32995 			 #  $use ExtraSamples
	DATATYPE = 32996 			 #  $use SampleFormat
	IMAGEDEPTH = 32997 			 #  z depth of image
	TILEDEPTH = 32998 			 #  z depth/data tile
# tags 33300-33309 are private tags registered to Pixar

# TIFFTAG_PIXAR_IMAGEFULLWIDTH and TIFFTAG_PIXAR_IMAGEFULLLENGTH
# are set when an image has been cropped out of a larger image.
# They reflect the size of the original uncropped image.
# The TIFFTAG_XPOSITION and TIFFTAG_YPOSITION can be used
# to determine the position of the smaller image in the larger one.

    PIXAR_IMAGEFULLWIDTH = 33300     # full image size in x
    PIXAR_IMAGEFULLLENGTH = 33301    # full image size in y

# tag 33405 is a private tag registered to Eastman Kodak
	WRITERSERIALNUMBER = 33405 			 #  device serial number
	CFAREPEATPATTERNDIM = 33421 			 #  dimensions of CFA pattern
	CFAPATTERN = 33422 			 #  color filter array pattern
# tag 33432 is listed in the 6.0 spec w/ unknown ownership
	COPYRIGHT = 33432 			 #  copyright string
# IPTC TAG from RichTIFF specifications
    RICHTIFFIPTC = 33723
# 34016-34029 are reserved for ANSI IT8 TIFF/IT <dkelly@apago.com)
	IT8SITE = 34016 			 #  site name
	IT8COLORSEQUENCE = 34017 			 #  color seq. [RGB,CMYK,etc]
	IT8HEADER = 34018 			 #  DDES Header
	IT8RASTERPADDING = 34019 			 #  raster scanline padding
	IT8BITSPERRUNLENGTH = 34020 			 #  # of bits in short run
	IT8BITSPEREXTENDEDRUNLENGTH = 34021 			 #  # of bits in long run
	IT8COLORTABLE = 34022 			 #  LW colortable
	IT8IMAGECOLORINDICATOR = 34023 			 #  BP/BL image color switch
	IT8BKGCOLORINDICATOR = 34024 			 #  BP/BL bg color switch
	IT8IMAGECOLORVALUE = 34025 			 #  BP/BL image color value
	IT8BKGCOLORVALUE = 34026 			 #  BP/BL bg color value
	IT8PIXELINTENSITYRANGE = 34027 			 #  MP pixel intensity value
	IT8TRANSPARENCYINDICATOR = 34028 			 #  HC transparency switch
	IT8COLORCHARACTERIZATION = 34029 			 #  color character. table
	IT8HCUSAGE = 34030 			 #  HC usage indicator
    IT8TRAPINDICATOR = 34031	 # Trapping indicator (untrapped=0, trapped=1)
	IT8CMYKEQUIVALENT = 34032 			 #  CMYK color equivalents
# tags 34232-34236 are private tags registered to Texas Instruments
	FRAMECOUNT = 34232 			 #  Sequence Frame Count
# tag 34377 is private tag registered to Adobe for PhotoShop
    PHOTOSHOP = 34377
# tags 34665, 34853 and 40965 are documented in EXIF specification
	EXIFIFD = 34665 			 #  Pointer to EXIF private directory
# tag 34750 is a private tag registered to Adobe?
	ICCPROFILE = 34675 			 #  ICC profile data
	IMAGELAYER = 34732 			 #  !TIFF/FX image layer information
# tag 34750 is a private tag registered to Pixel Magic
	JBIGOPTIONS = 34750 			 #  JBIG options
	GPSIFD = 34853 			 #  Pointer to GPS private directory
# tags 34908-34914 are private tags registered to SGI
	FAXRECVPARAMS = 34908 			 #  encoded Class 2 ses. parms
	FAXSUBADDRESS = 34909 			 #  received SubAddr string
	FAXRECVTIME = 34910 			 #  receive time (secs)
	FAXDCS = 34911 			 #  encoded fax ses. params, Table 2/T.30
# tags 37439-37443 are registered to SGI <gregl@sgi.com>
	STONITS = 37439 			 #  Sample value to Nits
# tag 34929 is a private tag registered to FedEx
    FEDEX_EDR = 34929	# unknown use
	INTEROPERABILITYIFD = 40965 			 #  Pointer to Interoperability private directory

	MODELPIXELSCALE = 33550
	MODELTIEPOINT = 33922

	GEOKEYDIRECTORY = 34735
	GEODOUBLEPARAMS = 34736
	GEOASCIIPARAMS = 34737

	GDALMETADATA = 42112
	GDALNODATA = 42113

# Adobe Digital Negative (DNG) format tags */
	DNGVERSION = 50706 			 #  &DNG version number
	DNGBACKWARDVERSION = 50707 			 #  &DNG compatibility version
	UNIQUECAMERAMODEL = 50708 			 #  &name for the camera model
	CFALAYOUT = 50711 			 #  &spatial layout of the CFA
	LINEARIZATIONTABLE = 50712 			 #  &lookup table description
	BLACKLEVEL = 50714 			 #  &zero light encoding level
	DEFAULTSCALE = 50718 			 #  &default scale factors
	CAMERACALIBRATION1 = 50723 			 #  &calibration matrix 1
	CAMERACALIBRATION2 = 50724 			 #  &calibration matrix 2
	BASELINENOISE = 50731 			 #  &relative noise level
	LINEARRESPONSELIMIT = 50734 			 #  &non-linear encoding range
	CAMERASERIALNUMBER = 50735 			 #  &camera's serial number
	LENSINFO = 50736 			 #  info about the lens
	CHROMABLURRADIUS = 50737 			 #  &chroma blur radius
#define TIFFTAG_ANTIALIASSTRENGTH	50738	&relative strength of the camera's anti-alias filter
	SHADOWSCALE = 50739 			 #  &used by Adobe Camera Raw
	DNGPRIVATEDATA = 50740 			 #  &manufacturer's private data
#define TIFFTAG_MAKERNOTESAFETY		50741	whether the EXIF MakerNote tag is safe to preserve along with the rest of the EXIF data
	CALIBRATIONILLUMINANT1 = 50778 			 #  &illuminant 1
	CALIBRATIONILLUMINANT2 = 50779 			 #  &illuminant 2
	BESTQUALITYSCALE = 50780 			 #  &best quality multiplier
#define TIFFTAG_RAWDATAUNIQUEID		50781	/* &unique identifier for
#define TIFFTAG_ORIGINALRAWFILENAME	50827	/* &file name of the original
#define TIFFTAG_ORIGINALRAWFILEDATA	50828	/* &contents of the original
#define TIFFTAG_ACTIVEAREA		50829	/* &active (non-masked) pixels
#define TIFFTAG_MASKEDAREAS		50830	/* &list of coordinates
	ASSHOTICCPROFILE = 50831 			 #  &these two tags used to
#define TIFFTAG_ASSHOTPREPROFILEMATRIX	50832	/* map cameras's color space
	CURRENTICCPROFILE = 50833 			 #  &
	CURRENTPREPROFILEMATRIX = 50834 			 #  &
	IMAGEJ_METADATA_BYTECOUNTS = 50838		# Embedded IMAGEJ metadata
	IMAGEJ_METADATA = 50839
	
# tag 65535 is an undefined tag used by Eastman Kodak
	DCSHUESHIFTVALUES = 65535 			 #  hue shift correction data

# The following are ``pseudo tags'' that can be used to control
# codec-specific functionality.  These tags are not written to file.
# Note that these values start at 0xffff+1 so that they'll never
# collide with Aldus-assigned tags.

# If you want your private pseudo tags ``registered'' (i.e. added to
# this file), please post a bug report via the tracking system at
# http://www.remotesensing.org/libtiff/bugs.html with the appropriate
# C definitions to add.

	FAXMODE = 65536 			 #  Group 3/4 format control
#define	    FAXMODE_CLASSIC	0x0000		/* default, include RTC */
#define	    FAXMODE_NORTC	0x0001		/* no RTC at end of data */
#define	    FAXMODE_NOEOL	0x0002		/* no EOL code at end of row */
#define	    FAXMODE_BYTEALIGN	0x0004		/* byte align row */
#define	    FAXMODE_WORDALIGN	0x0008		/* word align row */
#define	    FAXMODE_CLASSF	FAXMODE_NORTC	/* TIFF Class F */
	JPEGQUALITY = 65537 			 #  Compression quality level
# Note: quality level is on the IJG 0-100 scale.  Default value is 75
	JPEGCOLORMODE = 65538 			 #  Auto RGB<=>YCbCr convert?
#define	    JPEGCOLORMODE_RAW	0x0000		/* no conversion (default) */
#define	    JPEGCOLORMODE_RGB	0x0001		/* do auto conversion */
	JPEGTABLESMODE = 65539 			 #  What to put in JPEGTables
#define	    JPEGTABLESMODE_QUANT 0x0001		/* include quantization tbls */
#define	    JPEGTABLESMODE_HUFF	0x0002		/* include Huffman tbls */
# Note: default is JPEGTABLESMODE_QUANT | JPEGTABLESMODE_HUFF */
	FAXFILLFUNC = 65540 			 #  G3/G4 fill function
	PIXARLOGDATAFMT = 65549 			 #  PixarLogCodec I/O data sz
#define	    PIXARLOGDATAFMT_8BIT	0	/* regular u_char samples */
#define	    PIXARLOGDATAFMT_8BITABGR	1	/* ABGR-order u_chars */
#define	    PIXARLOGDATAFMT_11BITLOG	2	/* 11-bit log-encoded (raw) */
#define	    PIXARLOGDATAFMT_12BITPICIO	3	/* as per PICIO (1.0==2048) */
#define	    PIXARLOGDATAFMT_16BIT	4	/* signed short samples */
#define	    PIXARLOGDATAFMT_FLOAT	5	/* IEEE float samples */
# 65550-65556 are allocated to Oceana Matrix <dev@oceana.com> */
	DCSIMAGERTYPE = 65550 			 #  imager model & filter
#define     DCSIMAGERMODEL_M3           0       /* M3 chip (1280 x 1024) */
#define     DCSIMAGERMODEL_M5           1       /* M5 chip (1536 x 1024) */
#define     DCSIMAGERMODEL_M6           2       /* M6 chip (3072 x 2048) */
#define     DCSIMAGERFILTER_IR          0       /* infrared filter */
#define     DCSIMAGERFILTER_MONO        1       /* monochrome filter */
#define     DCSIMAGERFILTER_CFA         2       /* color filter array */
#define     DCSIMAGERFILTER_OTHER       3       /* other filter */
	DCSINTERPMODE = 65551 			 #  interpolation mode
#define     DCSINTERPMODE_NORMAL        0x0     /* whole image, default */
#define     DCSINTERPMODE_PREVIEW       0x1     /* preview of image (384x256) */
	DCSBALANCEARRAY = 65552 			 #  color balance values
	DCSCORRECTMATRIX = 65553 			 #  color correction values
	DCSGAMMA = 65554 			 #  gamma value
	DCSTOESHOULDERPTS = 65555 			 #  toe & shoulder points
	DCSCALIBRATIONFD = 65556 			 #  calibration file desc
# Note: quality level is on the ZLIB 1-9 scale. Default value is -1
	ZIPQUALITY = 65557 			 #  compression quality level
	PIXARLOGQUALITY = 65558 			 #  PixarLog uses same scale
# 65559 is allocated to Oceana Matrix <dev@oceana.com>
	DCSCLIPRECTANGLE = 65559 			 #  area of image to acquire
	SGILOGDATAFMT = 65560 			 #  SGILog user data format
#define     SGILOGDATAFMT_FLOAT		0	/* IEEE float samples */
#define     SGILOGDATAFMT_16BIT		1	/* 16-bit samples */
#define     SGILOGDATAFMT_RAW		2	/* uninterpreted data */
#define     SGILOGDATAFMT_8BIT		3	/* 8-bit RGB monitor values */
	SGILOGENCODE = 65561 			 #  SGILog data encoding control
#define     SGILOGENCODE_NODITHER	0     /* do not dither encoded values*/
#define     SGILOGENCODE_RANDITHER	1     /* randomly dither encd values */
	LZMAPRESET = 65562 			 #  LZMA2 preset (compression level)
	PERSAMPLE = 65563 			 #  interface for per sample tags
#define     PERSAMPLE_MERGED        0	/* present as a single value */
#define     PERSAMPLE_MULTI         1	/* present as multiple values */


# EXIF tags

	EXIFTAG_EXPOSURETIME = 33434 			 #  Exposure time
	EXIFTAG_FNUMBER = 33437 			 #  F number
	EXIFTAG_EXPOSUREPROGRAM = 34850 			 #  Exposure program
	EXIFTAG_SPECTRALSENSITIVITY = 34852 			 #  Spectral sensitivity
	EXIFTAG_ISOSPEEDRATINGS = 34855 			 #  ISO speed rating
	EXIFTAG_OECF = 34856 			 #  Optoelectric conversion factor
	EXIFTAG_EXIFVERSION = 36864 			 #  Exif version
	EXIFTAG_DATETIMEORIGINAL = 36867 			 #  Date and time of original data generation
	EXIFTAG_DATETIMEDIGITIZED = 36868 			 #  Date and time of digital data generation
	EXIFTAG_COMPONENTSCONFIGURATION = 37121 			 #  Meaning of each component
	EXIFTAG_COMPRESSEDBITSPERPIXEL = 37122 			 #  Image compression mode
	EXIFTAG_SHUTTERSPEEDVALUE = 37377 			 #  Shutter speed
	EXIFTAG_APERTUREVALUE = 37378 			 #  Aperture
	EXIFTAG_BRIGHTNESSVALUE = 37379 			 #  Brightness
	EXIFTAG_EXPOSUREBIASVALUE = 37380 			 #  Exposure bias
	EXIFTAG_MAXAPERTUREVALUE = 37381 			 #  Maximum lens aperture
	EXIFTAG_SUBJECTDISTANCE = 37382 			 #  Subject distance
	EXIFTAG_METERINGMODE = 37383 			 #  Metering mode
	EXIFTAG_LIGHTSOURCE = 37384 			 #  Light source
	EXIFTAG_FLASH = 37385 			 #  Flash
	EXIFTAG_FOCALLENGTH = 37386 			 #  Lens focal length
	EXIFTAG_SUBJECTAREA = 37396 			 #  Subject area
	EXIFTAG_MAKERNOTE = 37500 			 #  Manufacturer notes
	EXIFTAG_USERCOMMENT = 37510 			 #  User comments
	EXIFTAG_SUBSECTIME = 37520 			 #  DateTime subseconds
	EXIFTAG_SUBSECTIMEORIGINAL = 37521 			 #  DateTimeOriginal subseconds
	EXIFTAG_SUBSECTIMEDIGITIZED = 37522 			 #  DateTimeDigitized subseconds
	EXIFTAG_FLASHPIXVERSION = 40960 			 #  Supported Flashpix version
	EXIFTAG_COLORSPACE = 40961 			 #  Color space information
	EXIFTAG_PIXELXDIMENSION = 40962 			 #  Valid image width
	EXIFTAG_PIXELYDIMENSION = 40963 			 #  Valid image height
	EXIFTAG_RELATEDSOUNDFILE = 40964 			 #  Related audio file
	EXIFTAG_FLASHENERGY = 41483 			 #  Flash energy
	EXIFTAG_SPATIALFREQUENCYRESPONSE = 41484 			 #  Spatial frequency response
	EXIFTAG_FOCALPLANEXRESOLUTION = 41486 			 #  Focal plane X resolution
	EXIFTAG_FOCALPLANEYRESOLUTION = 41487 			 #  Focal plane Y resolution
	EXIFTAG_FOCALPLANERESOLUTIONUNIT = 41488 			 #  Focal plane resolution unit
	EXIFTAG_SUBJECTLOCATION = 41492 			 #  Subject location
	EXIFTAG_EXPOSUREINDEX = 41493 			 #  Exposure index
	EXIFTAG_SENSINGMETHOD = 41495 			 #  Sensing method
	EXIFTAG_FILESOURCE = 41728 			 #  File source
	EXIFTAG_SCENETYPE = 41729 			 #  Scene type
	EXIFTAG_CFAPATTERN = 41730 			 #  CFA pattern
	EXIFTAG_CUSTOMRENDERED = 41985 			 #  Custom image processing
	EXIFTAG_EXPOSUREMODE = 41986 			 #  Exposure mode
	EXIFTAG_WHITEBALANCE = 41987 			 #  White balance
	EXIFTAG_DIGITALZOOMRATIO = 41988 			 #  Digital zoom ratio
	EXIFTAG_FOCALLENGTHIN35MMFILM = 41989 			 #  Focal length in 35 mm film
	EXIFTAG_SCENECAPTURETYPE = 41990 			 #  Scene capture type
	EXIFTAG_GAINCONTROL = 41991 			 #  Gain control
	EXIFTAG_CONTRAST = 41992 			 #  Contrast
	EXIFTAG_SATURATION = 41993 			 #  Saturation
	EXIFTAG_SHARPNESS = 41994 			 #  Sharpness
	EXIFTAG_DEVICESETTINGDESCRIPTION = 41995 			 #  Device settings description
	EXIFTAG_SUBJECTDISTANCERANGE = 41996 			 #  Subject distance range
	EXIFTAG_IMAGEUNIQUEID = 42016 			 #  Unique image ID
end

@enum SampleFormats begin
	SAMPLEFORMAT_UINT = 1				# !unsigned integer data
	SAMPLEFORMAT_INT = 2				# !signed integer data
	SAMPLEFORMAT_IEEEFP = 3				# !IEEE floating point data 
	SAMPLEFORMAT_VOID = 4				# !untyped data
	SAMPLEFORMAT_COMPLEXINT = 5			# !complex signed int
	SAMPLEFORMAT_COMPLEXIEEEFP = 6		# !complex ieee floating
end

@enum PhotometricInterpretations begin
	PHOTOMETRIC_MINISWHITE = 0		# min value is white
	PHOTOMETRIC_MINISBLACK = 1		# min value is black
	PHOTOMETRIC_RGB = 2				# RGB color model
	PHOTOMETRIC_PALETTE = 3			# color map indexed
	PHOTOMETRIC_MASK = 4			# $holdout mask
	PHOTOMETRIC_SEPARATED = 5		# !color separations
	PHOTOMETRIC_YCBCR = 6			# !CCIR 601
	PHOTOMETRIC_CIELAB = 8			# !1976 CIE L*a*b*
	PHOTOMETRIC_ICCLAB = 9			# ICC L*a*b* [Adobe TIFF Technote 4]
	PHOTOMETRIC_ITULAB = 10			# ITU L*a*b*
	PHOTOMETRIC_CFA = 32803			# color filter array
	PHOTOMETRIC_LOGL = 32844		# CIE Log2(L)
	PHOTOMETRIC_LOGLUV = 32845		# CIE Log2(L) (u',v')
end

@enum CompressionType begin
	COMPRESSION_NONE = 1			# No compression
	COMPRESSION_CCITTRLE = 2		# CCITT modified Huffman RLE
	COMPRESSION_CCITT_T4 = 3        # CCITT T.4 (TIFF 6 name)
	COMPRESSION_CCITT_T6 = 4 		# CCITT T.6 (TIFF 6 name)
	COMPRESSION_LZW = 5				# Lempel-Ziv  & Welch
	COMPRESSION_OJPEG = 6			# !6.0 JPEG
	COMPRESSION_JPEG = 7			# %JPEG DCT compression
	COMPRESSION_T85	= 9				# !TIFF/FX T.85 JBIG compression
	COMPRESSION_T43 = 10			# !TIFF/FX T.43 colour by layered JBIG compression
	COMPRESSION_NEXT = 32766		# NeXT 2-bit RLE
	COMPRESSION_CCITTRLEW = 32771	#1 w/ word alignment
	COMPRESSION_PACKBITS = 32773	# Macintosh RLE
	COMPRESSION_THUNDERSCAN	= 32809	# ThunderScan RLE
# codes 32895-32898 are reserved for ANSI IT8 TIFF/IT <dkelly@apago.com)
	COMPRESSION_IT8CTPAD = 32895 	# IT8 CT w/padding
	COMPRESSION_IT8LW = 32896   	# IT8 Linework RLE
	COMPRESSION_IT8MP = 32897   	# IT8 Monochrome picture
	COMPRESSION_IT8BL = 32898   	# IT8 Binary line art
# compression codes 32908-32911 are reserved for Pixar
	COMPRESSION_PIXARFILM = 32908   # Pixar companded 10bit LZW
	COMPRESSION_PIXARLOG = 32909   	# Pixar companded 11bit ZIP
	COMPRESSION_DEFLATE = 32946		# Deflate compression
	COMPRESSION_ADOBE_DEFLATE = 8	# Deflate compression, as recognized by Adobe
# compression code 32947 is reserved for Oceana Matrix <dev@oceana.com>
	COMPRESSION_JBIG = 34661		# ISO JBIG
	COMPRESSION_SGILOG = 34676		# SGI Log Luminance RLE
	COMPRESSION_SGILOG24 = 34677	# SGI Log 24-bit packed
	COMPRESSION_JP2000 = 34712   	# Leadtools JPEG2000
	COMPRESSION_LZMA = 34925		# LZMA2
end

@enum ExtraSamples begin
	EXTRASAMPLE_UNSPECIFIED = 0			# unspecified data
	EXTRASAMPLE_ASSOCALPHA = 1			# associated alpha data
	EXTRASAMPLE_UNASSALPHA = 2			# unassociated alpha data
	# there's speculation online about this invalid code being used by CorelDraw
	EXTRASAMPLE_ASSOCALPHA_NS = 999		# associated alpha data (non-standard code)
end
