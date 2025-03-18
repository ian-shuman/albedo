PRO Albedo_Image_Statistics

;**********************************************************************************************;
; This function is building layer-stack of time-series data for albedo product


;           Latest updated on 2020/08/14 by Daryl Yang <<dediyang@bnl.gov>>
;**********************************************************************************************;
;Recrod the time of start
Start_Time=systime(1)

;*********************************** SET OUTPUT DIRECTORY *************************************;
out_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\map_analysis'
file_mkdir, out_dir
if out_dir then begin
  print, '...... Creating  "', out_dir, '" Successful'
endif else begin
  print, '...... Creating  "', out_dir, '" Failed'
endelse
;**********************************************************************************************;

;********************************** SET USER PARAMETER ****************************************;
window_size = 5
;**********************************************************************************************;

;*************************************** LODA DATA ********************************************;
; load transect roi
roi_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\Vegetation'
roi_file_dir = file_search(roi_dir, 'council_pft_map_reso30m_final.dat')
roi_raster = envi.OpenRaster(roi_file_dir[0])
roi_data = roi_raster.getData()

; load vegetation data
pft_cover_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\Vegetation\council_pft_fcover_reso30m_final.dat'
pft_cover_raster = envi.OpenRaster(pft_cover_dir)
; load topography data
topo_dem_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\Topography\54_16_2_1_2m_v3.0_reg_dem_utm_res30m_clipped.tif'
topo_dem_raster = envi.OpenRaster(topo_dem_dir)
topo_aspect_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\Topography\54_16_2_1_2m_v3.0_reg_aspect_utm_res30m_clipped.tif'
topo_aspect_raster = envi.OpenRaster(topo_aspect_dir)
topo_slope_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\Topography\54_16_2_1_2m_v3.0_reg_slope_utm_res30m_clipped.tif'
topo_slope_raster = envi.OpenRaster(topo_slope_dir)
topo_hillslope_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\Topography\54_16_2_1_2m_v3.0_reg_hillshade_utm_res30m_clipped.tif'
topo_hillslope_raster = envi.OpenRaster(topo_hillslope_dir)
; load albedo data
winter_albedo_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\albedo\winter_mean_albedo_v1.dat'
winter_albedo_raster = envi.OpenRaster(winter_albedo_dir)
summer_albedo_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\albedo\summer_mean_albedo_v1.dat'
summer_albedo_raster = envi.OpenRaster(summer_albedo_dir)
albedo_dir_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\albedo\winter_summer_albedo_difference_v1.dat'
albedo_dir_raster = envi.OpenRaster(albedo_dir_dir)
transition_date_dir = 'G:\MyWork\Albedo_Scaling\Data_Analysis\albedo\transition_date_v1.dat'
transition_date_raster = envi.OpenRaster(transition_date_dir)

data_combn_raster = EnviLayerStackRaster([roi_raster, pft_cover_raster, topo_dem_raster, topo_aspect_raster, topo_slope_raster, topo_hillslope_raster, $
  winter_albedo_raster, summer_albedo_raster, albedo_dir_raster, transition_date_raster])
;**********************************************************************************************;


n_column = roi_raster.ncolumns
n_row = roi_raster.nrows

transect_data_combn = []
for i=0, n_column-1 do begin
  print, i
  for j=0, n_row-1 do begin
    roi = roi_raster.GetData(sub_rect=[i,j,i,j])
    if roi gt 0 then begin
      ; extract vegetation data
      pixel_data = data_combn_raster.GetData(sub_rect=[i,j,i,j])
      pixel_data = reform(pixel_data)
      pixel_data_combn = [i,j, pixel_data]
      transect_data_combn = [[transect_data_combn],[pixel_data_combn]]
    endif
  endfor
endfor

headers = ['x', 'y', 'PFT', 'EVT', 'DTSA' , 'DTSW', 'DTSO', 'DLS', 'DDS', 'EVS', 'FOR', 'DGR', 'WGR', 'MOS',$
  'LIC', 'NVS', 'DEM', 'aspect' , 'slope', 'hillslope', 'winter_albedo', 'summer_albedo', 'albedo_dif', 'transition_date']

out_name = out_dir + '\transect_datav_v2.csv'
write_csv, out_name, transect_data_combn, header = headers
END