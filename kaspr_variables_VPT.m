% kaspr_variables_VPT.m
% Description: Read in KASPR VPT data.
% Author: Erin Leghart; erin.leghart@stonybrook.edu
% Last Updated: April 14, 2025

function [timeh, times, ref, spw, snr, rangekm, elev, verticalRes, file_duration_s] =...
    kaspr_variables_VPT(kasprdata)

    % Read the NetCDF File
    ncid = netcdf.open(kasprdata,'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid,'timeh'); timeh = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'times'); times = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'ref'); ref = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'spw'); spw = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'snr'); snr = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'rangekm'); rangekm = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'elev'); elev = netcdf.getVar(ncid,varid);
    verticalRes = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'range_resolution_m');
    verticalRes = str2double(verticalRes);
    file_duration_s = numel(timeh);
    netcdf.close(ncid);
end

