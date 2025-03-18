PRO STEP1_LANDSAT_CLIP
;;-------------------------------------------------------------------------------------------------------------;;
;;                         Lastest updated by Daryl YANG on Aug 26, 2018
;;   USAGE: this code is used for clipping images using a shape file
;;   INPUT PARAMETERS:
;;       path_IMG: the folder that contains the unzipped surface reflectance files
;;       vector_path: vector file directory
;;       path_OUT: the folder for storing the clipped files
;;-------------------------------------------------------------------------------------------------------------;;

; set the input parameters
path_IMG='F:\MyWork\Albedo_Scaling\Data_Processing\step1_Convert_HDF_to_GeoTIFF\2019'
file_list=file_search(path_IMG, '*.tif', count=num, /test_regular)

vector_path = '\\modex.bnl.gov\data2\dyang\projects\albedo_scaling\ngee_watersheds\boundary\teller.shp'
vectorFile = envi.OpenVector(vector_path)

path_OUT = '\\modex.bnl.gov\data2\dyang\projects\albedo_scaling\ngee_watersheds\teller\get_albedo\step2_clipped_albedo\2019'
FILE_MKDIR, path_OUT
if path_OUT then begin
  print, '--- Creating --- "', path_OUT, '" --- Successful ---'
endif else begin
  print, '--- Creating --- "', path_OUT, '" --- Failed ---'
endelse

; clip images
for i=0, num-1 do begin

  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN CONTINUE

  iIMG=file_list[i] 
  filename=getfilename(iIMG)
  filepathlen=strlen(filename)
  outname=STRMID(filename,0,filepathlen)
  
  ; create subarea based on vector file
  TASK = ENVITask('CreateSubrectsFromVector')
  rasterFile = envi.OpenRaster(iIMG)
  Task.INPUT_VECTOR = vectorFile
  Task.INPUT_RASTER = rasterFile
  Task.Execute
  
  ; pull out the subarea
  subArea = Task.SUBRECTS
  SubNames = Task.SUBRECT_NAMES
  
  ; clip the subarea
  DiceTask = ENVITask('DiceRasterBySubrects')
  DiceTask.INPUT_RASTER = rasterFile
  DiceTask.SUBRECT_ARRAY = subArea
  DiceTask.SUBRECT_NAMES = 'Council'
  DiceTask.Execute

  ; output the clipped images
  ExportTask = ENVITask('ExportRastersToDirectory')
  ExportTask.INPUT_RASTER = DiceTask.OUTPUT_RASTER
  ExportTask.OUTPUT_DIRECTORY = path_OUT
  ExportTask.Execute
endfor
END


Function getfilename, pathname
  ;The following functions are the subsidiary codes
  IF (N_PARAMS() NE 1) THEN RETURN, ''
  idx = STRPOS(pathname,'\', /REVERSE_SEARCH)
  filename = STRMID(pathname, idx + 1)
  RETURN, filename
END
