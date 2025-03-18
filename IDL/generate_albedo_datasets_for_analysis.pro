PRO Generate_Albedo_Datasets_for_Analysis

;**********************************************************************************************;
; This function is building layer-stack of time-series data for albedo product


;           Latest updated on 2020/08/14 by Daryl Yang <<dediyang@bnl.gov>>
;**********************************************************************************************;
;Recrod the time of start
Start_Time=systime(1)

;*********************************** SET OUTPUT DIRECTORY *************************************;
out_dir = '\\modex.bnl.gov\data2\dyang\projects\albedo_scaling\ngee_watersheds\teller\get_albedo\step5_seasonal_albedo'
file_mkdir, out_dir
if out_dir then begin
  print, '...... Creating  "', out_dir, '" Successful'
endif else begin
  print, '...... Creating  "', out_dir, '" Failed'
endelse
;**********************************************************************************************;

;********************************** SET USER PARAMETER ****************************************;
; day of year range for winter
winter_start = 1
winter_end = 121
; day of year range for summer
summer_start = 182
summer_end = 243
;**********************************************************************************************;

;*************************************** LODA DATA ********************************************;
; Load albedo file
data_dir = "\\modex.bnl.gov\data2\dyang\projects\albedo_scaling\ngee_watersheds\teller\get_albedo\step4_smoothed_ts"
file_dir = file_search(data_dir, "SG_Smoothed*")
; get raster data
data_raster = envi.OpenRaster(file_dir[0])
band_names = data_raster.metadata['band names']
albedo_data = data_raster.GetData()
; get map info
envi_open_file, file_dir[0], r_fid=fid
map_info=envi_get_map_info(fid=fid)
;**********************************************************************************************;

;********************************** CALCULATE SUMMER MEAN *************************************;
summer_albedo = albedo_data[*,*, where(band_names ge summer_start and band_names le summer_end)]
dims = size(summer_albedo, /dimensions)
; calculate mean summber albedo for each pixel
summer_albedo_mean = fltarr(dims[0], dims[1])
for i=0, dims[0]-1 do begin
  for j=0, dims[1]-1 do begin
    pixel_albedo = summer_albedo[i,j,*]
    ;loc = where(abs(pixel_albedo - mean(pixel_albedo, /NAN)) gt (stddev(pixel_albedo, /NAN)))
    ;pixel_albedo[loc] = 1 * !VALUES.F_NAN
    ;pixel_mean = mean(pixel_albedo, /NAN)
    resistant_mean, pixel_albedo, 1, pixel_mean
    summer_albedo_mean[i,j] = pixel_mean
  endfor
endfor
; save out summar albedo file
out_name = out_dir + '\' + 'summer_mean_albedo.dat'
envi_write_envi_file,summer_albedo_mean,out_name=out_name, map_info=map_info
;**********************************************************************************************;

;********************************** CALCULATE WINTER MEAN *************************************;
winter_albedo = albedo_data[*,*, where(band_names ge winter_start and band_names le winter_end)]
dims = size(winter_albedo, /dimensions)
; calculate mean summber albedo for each pixel
winter_albedo_mean = fltarr(dims[0], dims[1])
for i=0, dims[0]-1 do begin
  for j=0, dims[1]-1 do begin
    pixel_albedo = winter_albedo[i,j,*]
    ;loc = where(abs(pixel_albedo - mean(pixel_albedo, /NAN)) gt (stddev(pixel_albedo, /NAN)))
    ;pixel_albedo[loc] = 1 * !VALUES.F_NAN
    ;pixel_mean = mean(pixel_albedo, /NAN)
    resistant_mean, pixel_albedo, 1, pixel_mean
    winter_albedo_mean[i,j] = pixel_mean
  endfor
endfor
; save out summar albedo file
out_name = out_dir + '\' + 'winter_mean_albedo.dat'
envi_write_envi_file,winter_albedo_mean,out_name=out_name, map_info=map_info
;**********************************************************************************************;

;**************************** CALCULATE WINTER - SUMMER ***************************************;
winter_summer_dif = winter_albedo_mean-summer_albedo_mean
out_name = out_dir + '\' + 'winter_summer_albedo_diff.dat'
envi_write_envi_file,winter_summer_dif,out_name=out_name, map_info=map_info
;**********************************************************************************************;














END