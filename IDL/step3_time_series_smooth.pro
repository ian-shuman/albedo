PRO STEP3_TIME_SERIES_SMOOTH_V1
;**********************************************************************************************;
; This function is building layer-stack of time-series data for albedo product


;           Latest updated on 2020/08/14 by Daryl Yang <<dediyang@bnl.gov>>
;**********************************************************************************************;
;Recrod the time of start
Start_Time=systime(1)

;*********************************** SET OUTPUT DIRECTORY *************************************;
out_dir = '\\modex.bnl.gov\data2\dyang\projects\ngee_arctic\seward\analysis\phenology\council\planetscope\2022\VIs'
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
;**********************************************************************************************;

;*************************************** LODA DATA ********************************************;
; Load list of albedo image files
data_dir = '\\modex.bnl.gov\data2\dyang\projects\ngee_arctic\seward\analysis\phenology\council\planetscope\2022\VIs\snow_time_series.tif'
ts_raster = envi.OpenRaster(data_dir)
;band_names = ts_raster.metadata['band names']

date_stamp_orig = read_csv('\\modex.bnl.gov\data2\dyang\projects\ngee_arctic\seward\analysis\phenology\council\planetscope\2022\VIs\doy.csv', $
  header = header)
doy_stamp = date_stamp_orig.field2
band_names = doy_stamp
;**********************************************************************************************;
  
;********************************* SG filtering NDVI ******************************************;
ts_img = ts_raster.GetData()
ts_img = transpose(ts_img, [2, 1, 0])
ndvi_img_data = ts_img;/10000
dims_ndvi = size(ndvi_img_data, /dimensions)
ndvi_sg_filtered = fltarr(dims_ndvi)
for lon_col=0L, dims_ndvi[0]-1 do begin
  print, lon_col
  for lon_row = 0L, dims_ndvi[1]-1 do begin
    pixel_albedo_ts = reform(ndvi_img_data[lon_col, lon_row, *])
    pixel_albedo_ts_dropped = DropOutlier(pixel_albedo_ts, band_names)
    
    pixel_ndvi = smooth(pixel_albedo_ts_dropped, 15, /nan)
    
    pixel_ndvi[where(pixel_ndvi lt 0 or pixel_ndvi gt 1)] = 0
    pixel_ndvi[where(~finite(pixel_ndvi))] = 0
    subscribe_zero = where(pixel_ndvi gt 0 and pixel_ndvi lt 1, value_count)
    if (dims_ndvi[2] - value_count) lt dims_ndvi[2] then begin
      pixel_ndvi_filled = fill_up(temporary(pixel_ndvi), 1)
      ndvi_sg_filtered[lon_col, lon_row, *] = sgfilter(pixel_ndvi_filled)
    endif
  endfor
endfor

ndvi_sg_filtered_out = transpose(ndvi_sg_filtered, [1, 0, 2])
print, '...... saving SG smoothed NDVI time series'
envi_open_file, data_dir, r_fid=fid
map_info = envi_get_map_info(fid = fid)
file_name = file_basename(data_dir)
out_name = out_dir + '\' + 'SG_Smoothed_' + file_name
envi_write_envi_file, ndvi_sg_filtered_out, out_name = out_name, map_info = map_info, $
  bnames = ts_raster.metadata['BAND NAMES']
;**********************************************************************************************;  
print, 1
END


Function DropOutlier, ts_pixel, band_names
wd_size = 5
doy = fix(band_names[0:(n_elements(band_names)-2)])
doy_min = min(doy)
doy_max = max(doy)

doy_min2max = [doy_min:doy_max]

wd_doy_start = doy_min
wd_doy_end = wd_doy_start + wd_size

ts_pixel_updated = ts_pixel; fltarr(n_elements(ts_pixel))
while wd_doy_end le doy_max do begin
  wd_doy_loc = where(doy ge wd_doy_start and doy le wd_doy_end)
  wd_data = ts_pixel[wd_doy_loc]
  wd_mean = mean(wd_data, /NAN)
  wd_sd = stddev(wd_data, /NAN)
  wd_data[where((wd_data-wd_mean) gt wd_sd)] = 1 * !VALUES.F_NAN
  
  wd_doy_start = wd_doy_start + wd_size
  wd_doy_end = wd_doy_start + wd_size
  
  ts_pixel_updated[wd_doy_loc] = wd_data
  
endwhile
return, ts_pixel_updated
End






Function fill_up, vector_in, maxNDVI               ;remove the cloudy values

  num_elements=n_elements(vector_in)

  ;remove  continuous 0
  in_pixel = 0
  if vector_in[in_pixel] EQ 0 then begin
    pixel_start = in_pixel
    while vector_in[in_pixel] EQ 0 && in_pixel LT num_elements -2 do begin
      in_pixel = in_pixel +1
    endwhile
    pixel_end = in_pixel
    vector_in[pixel_start:pixel_end-1] = vector_in[pixel_end]
  endif

  in_pixel = num_elements-1
  if vector_in[in_pixel] EQ 0 then begin
    pixel_end = in_pixel
    while vector_in[in_pixel] EQ 0 && in_pixel GT 2 do begin
      in_pixel = in_pixel -1
    endwhile
    pixel_start = in_pixel
    vector_in[pixel_start+1:pixel_end] = vector_in[pixel_start]
  endif
  ;-----------
  in_pixel = 1
  while in_pixel LT num_elements - 2 do begin
    if vector_in[in_pixel] EQ 0 then begin
      pixel_start = in_pixel
      num_pixel = 1
      while vector_in[in_pixel] EQ 0 && in_pixel LT num_elements -2 do begin
        in_pixel = in_pixel +1
        num_pixel = num_pixel +1
      endwhile
      pixel_end = in_pixel
      temp = (vector_in[pixel_end] - vector_in[pixel_start-1]) / num_pixel
      for pixel = 0, num_pixel-2 do begin
        vector_in[pixel_start + pixel] = vector_in[pixel_start -1] + (pixel+1)*temp
      endfor
    endif
    in_pixel = in_pixel + 1
  endwhile

  in_pixel = 1
  while  in_pixel LT num_elements -2 do begin
    if (vector_in[in_pixel] - vector_in[in_pixel -1]) GE 0.2*maxNDVI $
      && (vector_in[in_pixel] - vector_in[in_pixel +1]) GE 0.2*maxNDVI then begin
      vector_in[in_pixel] = (vector_in[in_pixel -1] + vector_in[in_pixel +1]) /2.0
      in_pixel = in_pixel +2
    endif else begin
      in_pixel = in_pixel +1
    endelse
  endwhile

  return, vector_in

End

Function sgfilter,vector_in                 ;S-G filter

  num_elements = n_elements(vector_in)
  ; The first Savitzky-Golay fitting
  vector_in=reform(vector_in,num_elements)                         ; num_elements is the number of values of time-series
  savgolFilter = SAVGOL(15,15,0,2)                          ;set the window width(4,4) and degree (2) for computing trend curve
  rst = CONVOL(vector_in, savgolFilter, /EDGE_TRUNCATE)

  ; Calculate the threshold for loop control, so that the fit is maximize
  gdis = 0.0
  fl = IntARR(num_elements)

  for i =0,(num_elements-1) do begin
    fl[i] = (vector_in[i] ge rst[i])
    gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])*abs(vector_in[i]-rst[i])
  endfor

  ra4 = fltARR(num_elements)
  pre = fltARR(num_elements)

  ormax = gdis
  num   = 0

  loop_times = 0l
  while (gdis le ormax) && loop_times LT 10 do begin
    loop_times = loop_times +1
    for i =0,(num_elements-1) do begin
      ra4[i] = (vector_in[i] ge rst[i]) ? vector_in[i] : rst[i]
      pre[i] = rst[i]
    endfor

    ; The Savitzky-Golay fitting
    savgolFilter = SAVGOL(15, 15, 0, 3)        ;set the window width(4,4) and degree (6) for repetition
    rst = CONVOL(ra4, savgolFilter, /EDGE_TRUNCATE)
    ormax = gdis
    ; Calculate the fitting-effect index
    gdis = 0.0
    for i =0,(num_elements-1) do begin
      gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])*abs(vector_in[i]-rst[i])
    endfor
  endwhile

  if loop_times GE 1000 then begin
    print, 'loop times is: ', loop_times
  endif


  return, pre

End ; of function sgfilter













