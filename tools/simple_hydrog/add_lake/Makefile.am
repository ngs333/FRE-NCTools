#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
#
# FRE-NCTools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools.  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
bin_PROGRAMS = add_gages_from_table add_lake_data adj_gage_indices_reg \
	cp_small_lake_frac cr_big_lake_index_file cr_big_lakes  \
	cr_ground_type_files \
	cr_lake_files_gage cp_lake_files_new cr_lake_frac_files cr_lake_nc_from_idx \
	cr_lake_wbd_inputs get_basins  get_landf_diffs insert_lakes_3D insert_lakes \
	mod_ground_fracs mod_lake_indices

	##TODO:
	##cr_big_lakes_latlon

# Location of FMS library.
# TODO: allow builder to override and use a different libfms.
FMSROOT = $(top_builddir)/tools/simple_hydrog/libfmslite

AM_CFLAGS = $(NETCDF_CFLAGS)
AM_FCFLAGS = -I$(FMSROOT) $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_LIBS) $(NETCDF_FCLIBS)  $(RPATH_FLAGS)

add_gages_from_table_SOURCES = add_gages_from_table.f90
add_lake_data_SOURCES = add_lake_data.f90
adj_gage_indices_reg_SOURCES = adj_gage_indices_reg.f90
cp_small_lake_frac_SOURCES = cp_small_lake_frac.f90
cr_big_lake_index_file_SOURCES = cr_big_lake_index_file.f90
cr_big_lakes_SOURCES = cr_big_lakes.f90
#cr_big_lakes_latlon_SOURCES = cr_big_lakes_latlon.f90
cr_ground_type_files_SOURCES = cr_ground_type_files.f90
cr_lake_files_gage_SOURCES = cr_lake_files_gage.f90
#TODO cp_ vs cr_ for cp_lake_files_new - see the f90 file
cp_lake_files_new_SOURCES = cr_lake_files_new.f90
cr_lake_nc_from_idx_SOURCES = cr_lake_nc_from_idx.f90
cr_lake_wbd_inputs_SOURCES = cr_lake_wbd_inputs.f90
get_basins_SOURCES = get_basins.f90
cr_lake_frac_files_SOURCES = get_lake_cells.f90
get_landf_diffs_SOURCES = get_landf_diffs.f90
insert_lakes_3D_SOURCES = insert_lakes_3D.f90
insert_lakes_SOURCES = insert_lakes.f90
mod_ground_fracs_SOURCES = mod_ground_fracs.f90
mod_lake_indices_SOURCES = mod_lake_indices.f90

add_lake_data_LDADD =  $(FMSROOT)/libfms.a $(LDADD)
cp_small_lake_frac_LDADD =  $(FMSROOT)/libfms.a $(LDADD)
cr_big_lakes_LDADD =  $(FMSROOT)/libfms.a $(LDADD)
cr_lake_files_gage_LDADD =  $(FMSROOT)/libfms.a $(LDADD)
cp_lake_files_new_LDADD =  $(FMSROOT)/libfms.a $(LDADD)
cr_lake_nc_from_idx_LDADD =  $(FMSROOT)/libfms.a $(LDADD)
cr_lake_wbd_inputs_LDADD =  $(FMSROOT)/libfms.a $(LDADD)

# TODO: duplicate execution targets ?
#cr_ground_type_files_SOURCES = cr_ground_type_files_reyn.f90
#cr_lake_files_gage_SOURCES = cr_lake_files_c720.f90
#cr_lake_files_gage_SOURCES = get_lakes_large.f90
#cr_lake_nc_from_idx_SOURCES = cr_lake_nc_from_idx_chk_wbd.f90
#cr_lake_frac_files_SOURCES = cr_lake_frac_files_C48.f90


CLEANFILES = *.mod
