PRO STEP2_LAYERSTACK_IMAGES_V2
  ;**********************************************************************************************;
  ; This function is building layer-stack of time-series data for albedo product


  ;           Latest updated on 2020/08/14 by Daryl Yang <<dediyang@bnl.gov>>
  ;**********************************************************************************************;
  ;Recrod the time of start
  Start_Time=systime(1)

  ;*********************************** SET OUTPUT DIRECTORY *************************************;
  out_dir = '\\modex.bnl.gov\data2\dyang\projects\albedo_scaling\ngee_watersheds\council\get_albedo\step3_time_series'
  FILE_MKDIR, out_dir
  if out_dir then begin
    print, '...... Creating  "', out_dir, '" Successful'
  endif else begin
    print, '...... Creating  "', out_dir, '" Failed'
  endelse
  ;**********************************************************************************************;

  ;********************************** SET USER PARAMETER ****************************************;
  ; Day of year sequence
  date_use = indgen(1,366) + 1

  year = 2015
  ;**********************************************************************************************;

  ;*************************************** LODA DATA ********************************************;
  ; Load list of albedo image files
  data_dir = '\\modex.bnl.gov\data2\dyang\projects\albedo_scaling\ngee_watersheds\council\get_albedo\step2_clipped_albedo\' ; + strtrim(year, 1)
  img_list = file_search(data_dir, '*.dat')

  ;vector_path = 'G:\MyWork\Albedo_Scaling\ROIs\region_roi_shp.shp'
  ;vectorFile = envi.OpenVector(vector_path)

  ; Load benchmark image file to determine the spatial extent
  ref_raster = '\\modex.bnl.gov\data2\dyang\projects\albedo_scaling\ngee_watersheds\council\get_albedo\ref_raster\council_ref.dat'
  raster0 = envi.OpenRaster(ref_raster)
  ;**********************************************************************************************;

  ;********************************* EXTRACT FILE YY:MM:DD  *************************************;
  acquisition_date = strarr(2, n_elements(img_list))
  for i=0, n_elements(img_list)-1 do begin
    file = img_list[i]
    filename = getfilename(file)
    pattern = STRMID(filename, 19, 6)
    idoy = date_to_doy(pattern, doy, yr)
    acquisition_date[*, i] = [pattern, doy]
  endfor
  ;**********************************************************************************************;

  ;******************************** LAYSTACK FILES IN ORDER  ************************************;
  num_dates = n_elements(date_use)
  layerStack = raster0
  band_names = []
  for i=0,  num_dates-1 do begin ;num_dates-1
    print, i
    dat = date_use[i]
    matched_loc = where(acquisition_date[1,*] eq dat)
    ;WHERE(STRMATCH(acquisition_date[1], '*' + dat + '*', /FOLD_CASE) EQ 1)
    doy_layerStack = raster0
    if matched_loc[0] ne -1 then begin
      matched_file = img_list[matched_loc]
      for j=0, n_elements(matched_file)-1 do begin

        catch, error_status
        ;;This statement begins the error handler:
        if error_status ne 0 then continue
        print, 'error', error_status

        file_dir = matched_file[j]
        file_raster = envi.OpenRaster(file_dir)
        doy_layerStack = ENVILayerStackRaster([doy_layerStack, file_raster])
        file_raster.close
      endfor
      doy_rasters = doy_layerStack.GetData(bands=[1:(doy_layerStack.NBANDS-1)])
      
      doy_raste_update = intarr(size(doy_rasters, /dimension)) * !VALUES.F_NAN
      for j=0, doy_layerStack.NBANDS-2 do begin
        mask = intarr(doy_layerStack.ncolumns, doy_layerStack.nrows) * !VALUES.F_NAN
        band = doy_rasters[*,*,j]
        mask[where(band gt 1 and band le 15000)] = 1
        TEMP = band * mask
        doy_raste_update[*,*,j] = TEMP
      endfor

      temp = doy_layerStack.NBANDS-1
      if temp gt 1 then begin
        doy_mean = mean(doy_raste_update, dimension = 3, /NAN)
      endif
      if temp eq 1 then begin
        doy_mean = doy_raste_update
      endif
      
      doy_newRaster = ENVIRaster(doy_mean, URI=newFile, SPATIALREF = doy_layerStack.spatialref)
      doy_newRaster.Save
      
      layerStack = ENVILayerStackRaster([layerStack, doy_newRaster])
      band_names = [band_names, strtrim(dat)]
      
      doy_newRaster.close
    endif
  endfor

  data_combn = layerStack.GetData(BANDS=[1:(layerStack.NBANDS-1)])
  ; Replace nan values

  ;for i=0, layerStack.NBANDS-2 do begin
  ;  band = data_combn[*,*,i]
 ;   band[where(band gt 15000)] = 1*!VALUES.F_NAN
 ;   data_combn[*,*,i] = band
 ; endfor

  out_name = out_dir + '\' + 'council_lidar_area_albedo_ts_2015_2019.dat'; + strtrim(year, 1) + '.dat'
  newRaster = ENVIRaster(data_combn, URI=out_name, SPATIALREF = layerStack.spatialref)
  metadata = newRaster.METADATA
  metadata['band names'] = band_names
  newRaster.Save

END




Function getfilename, pathname
  ;The following functions are the subsidiary codes
  IF (N_PARAMS() NE 1) THEN RETURN, ''
  idx = STRPOS(pathname,'\', /REVERSE_SEARCH)
  filename = STRMID(pathname, idx + 1)
  RETURN, filename
End

Function date_to_doy, idate, DOY, yr
  ; ------------------------------------------------------------
  ;+              18-Sep-91
  ; NAME:
  ;   Date2DOY
  ; PURPOSE:
  ;   Convert yymmdd into DOY (day of the year).  Input can
  ;   be either string or integer.  The year is an Optional
  ;   return parameter.
  ; CALLING SEQUENCE:
  ;   Date2DOY, idate, DOY [, yr]
  ; INPUT:
  ;   idate input format for the date: yymmdd.
  ;     Data-type can be either string or integer.
  ; OUTPUT:
  ;   DOY integer with the day of the year.
  ; Output/Optional:
  ;   yr  year of the returned DOY.
  ; Note: If input data-type is string the returned values are
  ;   string-type and if input-type is longword the returned
  ;   parameters (DOY and yr) are integers.
  ; HISTORY:
  ;   written by GAL 18-Sep-91
  ;-
  ; -----------------------------------------------------------------
  ; ON_ERROR, 2 ;force a return to caller on error

  ; Check data type of input set ascII flag and convert to yy,mm,dd:
  info = SIZE(idate)
  IF (info(0) eq 0) THEN BEGIN
    scalar = 1        ;scalar flag set
  ENDIF ELSE BEGIN
    scalar = 0        ;vector input
  ENDELSE

  IF (info(info(0) + 1) eq 7) THEN BEGIN
    ascII = 1       ;ascII input flag set
    yy = FIX(STRMID(idate,0,2))   ;extract year
    mm = FIX(STRMID(idate,2,2))   ;extract month
    dd = FIX(STRMID(idate,4,2))   ;extract day
  ENDIF ELSE BEGIN      ;should be a longWord
    ascII = 0       ;non-ascII input
    sdate = STRTRIM(STRING(idate),2)  ;convert to string
    yy = FIX(STRMID(sdate,0,2))   ;extract year
    mm = FIX(STRMID(sdate,2,2))   ;extract month
    dd = FIX(STRMID(sdate,4,2))   ;extract day
  ENDELSE

  ; Check for leap year and compute DOY:
  ;               J   F   M   A   M   J   J   A   S   O   N   D
  imonth = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

  IF (scalar) THEN BEGIN      ;scalar input
    IF ((yy MOD 4) eq 0) THEN BEGIN ;leap year
      imonth(2) = 29      ;set feb
    ENDIF
    DOY = FIX( TOTAL(imonth(0:mm-1)) ) + dd
  ENDIF ELSE BEGIN
    DOY = dd        ;set correct len on vector
    leapYrs = WHERE( (yy MOD 4) eq 0) ;index of leap years
    nonLeap = WHERE( (yy MOD 4) ne 0) ;index of non-leap years
    IF (nonLeap(0) ne -1) THEN BEGIN
      FOR i=0, N_elements(nonLeap)-1 DO BEGIN
        DOY[nonLeap(i)] = FIX(TOTAL(imonth(0:mm(nonLeap(i))-1))) + dd(nonLeap(i))
      ENDFOR
    ENDIF
    IF (leapYrs(0) ne -1) THEN BEGIN
      imonth(2) = 29      ;set feb
      FOR i =0, N_elements(leapYrs)-1 DO BEGIN
        DOY[leapYrs(i)] = FIX(TOTAL(imonth(0:mm(leapYrs(i))-1))) + dd(leapYrs(i))
      ENDFOR
    ENDIF
  ENDELSE

  IF (N_PARAMS() EQ 3) THEN BEGIN         ;pass year back to caller
    IF (ascII) THEN BEGIN
      DOY = STRTRIM( STRING(DOY), 2)  ;convert to string
      yr = STRTRIM( STRING(yy), 2)  ;convert to string
    ENDIF ELSE BEGIN
      yr = yy
    ENDELSE
  ENDIF ELSE BEGIN      ;pass DOY only
    IF (ascII) THEN BEGIN
      DOY = STRTRIM( STRING(DOY), 2)  ;convert to string
    ENDIF
  ENDELSE

End

