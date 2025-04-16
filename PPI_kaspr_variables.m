% PPI_kaspr_variables.m
% Description: Read in KASPR PPI data.
% Author: Erin Leghart; erin.leghart@stonybrook.edu
% Last Updated: April 14, 2025

function [timeh, times, ref, spw, snr, rangekm, xkm, ykm, zkm, elev_deg, az_deg, file_duration_s] =...
    PPI_kaspr_variables(kasprdata)

    % Read the NetCDF File
    ncid = netcdf.open(kasprdata,'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid,'timeh'); timeh = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'times'); times = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'ref'); ref = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'spw'); spw = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'snr'); snr = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'rangekm'); rangekm = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'xkm'); xkm = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'ykm'); ykm = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'zkm'); zkm = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'elev_deg'); elev_deg = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'az_deg'); az_deg = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'file_duration_s'); file_duration_s = netcdf.getVar(ncid,varid);
    netcdf.close(ncid);
end
