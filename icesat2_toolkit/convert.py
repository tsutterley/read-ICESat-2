"""
convert.py
Written by Tyler Sutterley (12/2022)
Utilities for converting ICESat-2 HDF5 files into different formats

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org
        http://docs.h5py.org/en/stable/index.html
    zarr: Chunked, compressed, N-dimensional arrays in Python
        https://github.com/zarr-developers/zarr-python
        https://zarr.readthedocs.io/en/stable/index.html
    pandas: Python Data Analysis Library
        https://pandas.pydata.org/

UPDATE HISTORY:
    Updated 12/2022: place some imports behind try/except statements
    Updated 06/2022: place zarr and pandas imports behind try/except statements
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 01/2022: added ascii and dataframe outputs for ATL07
    Updated 09/2021: added ground track and time to output dataframes
        calculate a global reference point for ATL06/07/08 dataframes
    Updated 07/2021: comment for column number in yaml headers
    Updated 10/2020: added ascii output for ATL08
    Updated 08/2020: added output in pandas dataframe for ATL06 and ATL08
    Written 06/2020
"""
import re
import pathlib
import warnings
import itertools
import posixpath
import numpy as np

# attempt imports
try:
    import h5py
except ModuleNotFoundError:
    warnings.warn("h5py not available", ImportWarning)
try:
    import pandas
except ModuleNotFoundError:
    warnings.warn("pandas not available", ImportWarning)
try:
    import zarr
except ModuleNotFoundError:
    warnings.warn("zarr not available", ImportWarning)

class convert():
    np.seterr(invalid='ignore')
    def __init__(self, filename=None, reformat=None):
        """Utilities for converting ICESat-2 HDF5 files into different formats

        Parameters
        ----------
        filename: str, obj or NoneType, default None
            input HDF5 filename or io.BytesIO object
        reformat: str or NoneType
            output format

                - ``'csv'``: comma-separated values for each beam group
                - ``'dataframe'``: Pandas dataframe
                - ``'HDF5'``: rechunked HDF5
                - ``'txt'``: tab-delimited ascii for each beam group
                - ``'zarr'``: chunked zarr
        """
        self.filename = filename
        self.reformat = reformat
        # gps-based epoch for delta times #
        self.atlas_sdp_epoch = np.datetime64('2018-01-01T00:00:00')

    # PURPOSE: wrapper function for converting HDF5 files to another type
    def file_converter(self, **kwds):
        """
        Convert a HDF5 file to another format

        Parameters
        ----------
        **kwds: dict
            keyword arguments for output converter
        """
        if (self.reformat == 'zarr'):
            # output zarr file
            self.HDF5_to_zarr(**kwds)
        elif (self.reformat == 'HDF5'):
            # output rechunked HDF5 file
            self.HDF5_to_HDF5(**kwds)
        # elif (reformat == 'JPL'):
        #     # output JPL captoolkit formatted HDF5 files
        #     self.HDF5_to_JPL_HDF5(**kwds)
        elif self.reformat in ('csv','txt'):
            # output reduced files to ascii formats
            self.HDF5_to_ascii(**kwds)
        elif self.reformat in ('dataframe',):
            # output reduced files to pandas dataframe
            return self.HDF5_to_dataframe(**kwds)
        else:
            raise ValueError(f'Unknown format {self.reformat}')

    # PURPOSE: convert the HDF5 file to zarr copying all file data
    def HDF5_to_zarr(self, **kwds):
        """
        convert a HDF5 file to zarr copying all file data

        Parameters
        ----------
        **kwds: dict
            keyword arguments for output zarr converter
        """
        # output zarr file
        if isinstance(self.filename, (str, pathlib.Path)):
            zarr_file = pathlib.Path(self.filename).with_suffix('.zarr')
        else:
            zarr_file = pathlib.Path(self.filename.filename).with_suffix('.zarr')
        # copy everything from the HDF5 file to the zarr file
        with h5py.File(self.filename, mode='r') as source:
            dest = zarr.open_group(zarr_file, mode='w')
            # value checks on output zarr
            if not hasattr(dest, 'create_dataset'):
                raise ValueError('dest must be a group, got {!r}'.format(dest))
            # for each key in the root of the hdf5 file structure
            for k in source.keys():
                self.copy_from_HDF5(source[k], dest, name=k, **kwds)

    # PURPOSE: rechunk the HDF5 file copying all file data
    def HDF5_to_HDF5(self, **kwds):
        """
        rechunk a HDF5 file copying all file data

        Parameters
        ----------
        **kwds: dict
            keyword arguments for output HDF5 converter
        """
        # output HDF5 file
        if isinstance(self.filename, (str, pathlib.Path)):
            hdf5_file = pathlib.Path(self.filename).with_suffix('.h5')
        else:
            hdf5_file = pathlib.Path(self.filename.filename).with_suffix('.h5')
        # copy everything from the HDF5 file
        with h5py.File(self.filename,mode='r') as source:
            dest = h5py.File(hdf5_file, mode='w')
            # value checks on output HDF5
            if not hasattr(dest, 'create_dataset'):
                raise ValueError('dest must be a group, got {!r}'.format(dest))
            # for each key in the root of the hdf5 file structure
            for k in source.keys():
                self.copy_from_HDF5(source[k], dest, name=k, **kwds)

    # PURPOSE: Copy a named variable from the HDF5 file to the destination file
    def copy_from_HDF5(self, source, dest, name=None, **kwds):
        """
        Copy a named variable from the ``source`` HDF5 into the ``dest`` file

        Parameters
        ----------
        source: obj
            open file object for input
        dest: obj
            open file object for output
        name: str or NoneType, default None
            variable or group name
        chunks: int
            chunk size for output
        """
        if hasattr(source, 'shape') and bool(source.chunks):
            # if data can be chunked
            # copy a dataset/array
            if dest is not None and name in dest:
                raise RuntimeError('an object {!r} already exists in '
                    'destination {!r}'.format(name, dest.name))
            # setup creation keyword arguments
            create_kwds = {k:v for k,v in kwds.items() if v}
            if 'chunks' in create_kwds.keys():
                # setup chunks option to limit by dimensions
                chunks = ()
                for d in range(source.ndim):
                    chunks += (min([create_kwds.get('chunks'),source.shape[d]]),)
                create_kwds['chunks'] = chunks
            else:
                # setup chunks option, preserve by default
                chunks = ()
                for d in range(source.ndim):
                    chunks += (min([source.chunks[d],source.shape[d]]),)
                create_kwds.setdefault('chunks', chunks)
            # setup compression options
            # from h5py to zarr: use zarr default compression options
            if (self.reformat == 'zarr') and source.fillvalue:
                create_kwds.setdefault('fill_value', source.fillvalue)
            elif (self.reformat == 'HDF5') and source.fillvalue:
                create_kwds.setdefault('fillvalue', source.fillvalue)
            # create new dataset in destination
            ds = dest.create_dataset(name, shape=source.shape,
                dtype=source.dtype, **create_kwds)
            # copy data going chunk by chunk to avoid loading in entirety
            shape = ds.shape
            chunks = ds.chunks
            chunk_offsets = [range(0, s, c) for s, c in zip(shape, chunks)]
            for offset in itertools.product(*chunk_offsets):
                sel = tuple(slice(o, min(s, o + c)) for o, s, c in
                    zip(offset, shape, chunks))
                ds[sel] = source[sel]
            # copy attributes
            attrs = {key:self.attributes_encoder(source.attrs[key]) for key in
                source.attrs.keys() if self.attributes_encoder(source.attrs[key])}
            # update chunks attribute
            attrs['_ChunkSizes'] = self.attributes_encoder(create_kwds['chunks'])
            ds.attrs.update(attrs)
        elif hasattr(source, 'shape'):
            # if data cannot be chunked
            # copy a dataset/array
            if dest is not None and name in dest:
                raise RuntimeError('an object {!r} already exists in '
                    'destination {!r}'.format(name, dest.name))
            # setup creation keyword arguments
            create_kwds = {}
            # setup compression options
            # from h5py to zarr: use zarr default compression options
            if (self.reformat == 'zarr') and source.fillvalue:
                create_kwds.setdefault('fill_value', source.fillvalue)
            elif (self.reformat == 'HDF5') and source.fillvalue:
                create_kwds.setdefault('fillvalue', source.fillvalue)
            # create new dataset in destination
            ds = dest.create_dataset(name, shape=source.shape,
                dtype=source.dtype, **create_kwds)
            ds[:] = source[:]
            # copy attributes
            attrs = {key:self.attributes_encoder(source.attrs[key]) for key in
                source.attrs.keys() if self.attributes_encoder(source.attrs[key])}
            ds.attrs.update(attrs)
        else:
            # copy a group
            if (dest is not None and name in dest and hasattr(dest[name], 'shape')):
                raise RuntimeError('an array {!r} already exists in '
                    'destination {!r}'.format(name, dest.name))
            # require group in destination
            grp = dest.require_group(name)
            # copy attributes
            attrs = {key:self.attributes_encoder(source.attrs[key]) for key in
                source.attrs.keys() if self.attributes_encoder(source.attrs[key])}
            grp.attrs.update(attrs)
            # recursively copy from source
            for k in source.keys():
                self.copy_from_HDF5(source[k], grp, name=k, **kwds)

    # PURPOSE: reduce HDF5 files to beam groups and output to ascii
    def HDF5_to_ascii(self, **kwds):
        """
        Convert a HDF5 file to beam-level ascii files copying reduced sets of data

        Parameters
        ----------
        **kwds: dict
            keyword arguments for output ascii converter
        """
        # compile regular expression operator for extracting info from ICESat2 files
        rx = re.compile(r'(processed_)?(ATL\d+)(-\d{2})?_(\d{4})(\d{2})(\d{2})'
            r'(\d{2})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
        # split extension from HDF5 file
        # extract parameters from ICESat2 HDF5 file
        if isinstance(self.filename, (str, pathlib.Path)):
            hdf5_file = pathlib.Path(self.filename)
        else:
            hdf5_file = pathlib.Path(self.filename.filename)
        # extract parameters from ICESat2 HDF5 file
        SUB,PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX = \
            rx.findall(hdf5_file.name).pop()
        # output file suffix for csv or tab-delimited text
        delimiter = ',' if self.reformat == 'csv' else '\t'
        # copy bare minimum variables from the HDF5 file to the ascii file
        source = h5py.File(self.filename, mode='r')

        # find valid beam groups by testing for particular variables
        if (PRD == 'ATL06'):
            VARIABLE_PATH = ['land_ice_segments','segment_id']
        elif (PRD == 'ATL07'):
            VARIABLE_PATH = ['sea_ice_segments','height_segment_id']
        elif (PRD == 'ATL08'):
            VARIABLE_PATH = ['land_segments','segment_id_beg']
        elif (PRD == 'ATL10'):
            VARIABLE_PATH = ['freeboard_beam_segments','delta_time']
        elif (PRD == 'ATL12'):
            VARIABLE_PATH = ['ssh_segments','delta_time']
        # create list of valid beams within the HDF5 file
        beams = []
        for gtx in [k for k in source.keys() if bool(re.match(r'gt\d[lr]',k))]:
            # check if subsetted beam contains data
            try:
                source['/'.join([gtx,*VARIABLE_PATH])]
            except KeyError:
                pass
            else:
                beams.append(gtx)

        # for each valid beam within the HDF5 file
        for gtx in sorted(beams):
            # extract variables and attributes for each variable
            values = {}
            vattrs = {}
            # create a column stack of valid output segment values
            if (PRD == 'ATL06'):
                # land ice height
                var = source[gtx]['land_ice_segments']
                valid, = np.nonzero(var['atl06_quality_summary'][:] == 0)
                # variables for the output ascii file
                vnames = ['segment_id','delta_time','latitude','longitude',
                    'h_li','h_li_sigma']
                vformat = ('{1:0.0f}{0}{2:0.9f}{0}{3:0.9f}{0}{4:0.9f}{0}'
                    '{5:0.9f}{0}{6:0.9f}')
                # for each output variable
                for i,v in enumerate(vnames):
                    # convert data to numpy array for HDF5 compatibility
                    values[v] = np.copy(var[v][:])
                    # extract attributes
                    vattrs[v] = {atn:atv for atn,atv in var[v].attrs.items()}
                    # add precision attributes for ascii yaml header
                    if (v == 'segment_id'):
                        vattrs[v]['precision'] = 'integer'
                        vattrs[v]['units'] = 'count'
                    else:
                        vattrs[v]['precision'] = 'double_precision'
                    vattrs[v]['comment'] = f'column {i+1:d}'
            elif (PRD == 'ATL07'):
                # sea ice height
                var = source[gtx]['sea_ice_segments']
                valid, = np.nonzero(var['heights/height_segment_quality'][:] == 1)
                # variables for the output ascii file
                vnames = ['height_segment_id','delta_time',
                    'latitude','longitude','seg_dist_x',
                    'heights/height_segment_height',
                    'heights/height_segment_confidence',
                    'heights/height_segment_type',
                    'heights/height_segment_ssh_flag',
                    'heights/height_segment_w_gaussian',
                    'stats/photon_rate','stats/cloud_flag_asr',
                    'geophysical/height_segment_lpe',
                    'geophysical/height_segment_mss',
                    'geophysical/height_segment_ocean',
                    'geophysical/height_segment_ib']
                vformat = ('{1:0.0f}{0}{2:0.9f}{0}{3:0.9f}{0}{4:0.9f}{0}'
                    '{5:0.9f}{0}{6:0.9f}{0}{7:0.9f}{0}{8:0.0f}{0}{9:0.0f}{0}'
                    '{10:0.9f}{0}{11:0.9f}{0}{12:0.0f}{0}{13:0.9f}{0}'
                    '{14:0.9f}{0}{15:0.9f}{0}{16:0.9f}')
                # for each output variable
                for i,v in enumerate(vnames):
                    # convert data to numpy array for HDF5 compatibility
                    values[v] = np.copy(var[v][:])
                    # extract attributes
                    vattrs[v] = {atn:atv for atn,atv in var[v].attrs.items()}
                    # add precision attributes for ascii yaml header
                    if v in ('height_segment_id','heights/height_segment_type',
                             'heights/height_segment_ssh_flag',
                             'stats/cloud_flag_asr'):
                        vattrs[v]['precision'] = 'integer'
                    else:
                        vattrs[v]['precision'] = 'double_precision'
                    vattrs[v]['comment'] = f'column {i+1:d}'
            elif (PRD == 'ATL08'):
                # land and vegetation height
                var = source[gtx]['land_segments']
                valid, = np.nonzero(var['terrain/h_te_best_fit'][:] !=
                    var['terrain/h_te_best_fit'].fillvalue)
                # variables for the output ascii file
                vnames = ['segment_id_beg','segment_id_end','delta_time',
                    'latitude','longitude','terrain/h_te_best_fit',
                    'terrain/h_te_uncertainty','terrain/terrain_slope',
                    'canopy/h_canopy','canopy/h_canopy_uncertainty']
                vformat = ('{1:0.0f}{0}{2:0.0f}{0}{3:0.9f}{0}{4:0.9f}{0}'
                    '{5:0.9f}{0}{6:0.9f}{0}{7:0.9f}{0}{8:0.9f}{0}{9:0.9f}{0}'
                    '{10:0.9f}')
                # for each output variable
                for i,v in enumerate(vnames):
                    # convert data to numpy array for HDF5 compatibility
                    values[v] = np.copy(var[v][:])
                    # extract attributes
                    vattrs[v] = {atn:atv for atn,atv in var[v].attrs.items()}
                    # add precision attributes for ascii yaml header
                    if v in ('segment_id_beg','segment_id_end'):
                        vattrs[v]['precision'] = 'integer'
                        vattrs[v]['units'] = 'count'
                    else:
                        vattrs[v]['precision'] = 'double_precision'
                    vattrs[v]['comment'] = f'column {i+1:d}'

            # column stack of valid output segment values
            output = np.column_stack([values[v][valid] for v in vnames])

            # output ascii file
            granule = f'{hdf5_file.stem}_{gtx}.{self.reformat}'
            ascii_file = hdf5_file.parent.joinpath(granule)
            fid = ascii_file.open(mode='w', encoding='utf8')
            # print YAML header to top of file
            fid.write('{0}:\n'.format('header'))
            # global attributes for file
            fid.write('  {0}:\n'.format('global_attributes'))
            for atn,atv in source.attrs.items():
                if atn not in ('Conventions','Processing Parameters','hdfversion',
                    'history','identifier_file_uuid'):
                    fid.write('    {0:22}: {1}\n'.format(atn,self.attributes_encoder(atv)))
            # beam attributes
            fid.write('\n  {0}:\n'.format('beam_attributes'))
            for atn,atv in source[gtx].attrs.items():
                if atn not in ('Description',):
                    fid.write('    {0:22}: {1}\n'.format(atn,self.attributes_encoder(atv)))
            # data dimensions
            fid.write('\n  {0}:\n'.format('dimensions'))
            nrow,ncol = np.shape(output)
            fid.write('    {0:22}: {1:d}\n'.format('segments',nrow))
            # non-standard attributes
            fid.write('\n  {0}:\n'.format('non-standard_attributes'))
            # value to convert to GPS seconds (seconds since 1980-01-06T00:00:00)
            fid.write('    {0:22}:\n'.format('atlas_sdp_gps_epoch'))
            atlas_sdp_gps_epoch = source['ancillary_data']['atlas_sdp_gps_epoch']
            for atn in ['units','long_name']:
                atv = self.attributes_encoder(atlas_sdp_gps_epoch.attrs[atn])
                fid.write('      {0:20}: {1}\n'.format(atn,atv))
            fid.write('      {0:20}: {1:0.0f}\n'.format('value',atlas_sdp_gps_epoch[0]))
            # print variable descriptions to YAML header
            fid.write('\n  {0}:\n'.format('variables'))
            for v in vnames:
                fid.write('    {0:22}:\n'.format(posixpath.basename(v)))
                for atn in ['precision','units','long_name','comment']:
                    atv = self.attributes_encoder(vattrs[v][atn])
                    fid.write('      {0:20}: {1}\n'.format(atn,atv))
            # end of header
            fid.write('\n\n# End of YAML header\n')
            # print data to file
            for row in output:
                print(vformat.format(delimiter,*row),file=fid)
            # close the file
            fid.close()
        # close the source HDF5 file
        source.close()

    # PURPOSE: reduce HDF5 files to pandas dataframe
    def HDF5_to_dataframe(self, **kwds):
        """
        Convert a HDF5 file to a pandas dataframe copying reduced sets of data

        Parameters
        ----------
        **kwds: dict
            keyword arguments for output dataframe converter
        """
        # compile regular expression operator for extracting info from ICESat2 files
        rx = re.compile(r'(processed_)?(ATL\d+)(-\d{2})?_(\d{4})(\d{2})(\d{2})'
            r'(\d{2})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
        # split extension from HDF5 file
        # extract parameters from ICESat2 HDF5 file
        if isinstance(self.filename, (str, pathlib.Path)):
            hdf5_file = pathlib.Path(self.filename)
        else:
            hdf5_file = pathlib.Path(self.filename.filename)
        # extract parameters from ICESat2 HDF5 file
        SUB,PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX = \
            rx.findall(hdf5_file.name).pop()

        # copy bare minimum variables from the HDF5 file to pandas data frame
        source = h5py.File(self.filename,mode='r')

        # find valid beam groups by testing for particular variables
        if (PRD == 'ATL06'):
            VARIABLE_PATH = ['land_ice_segments','segment_id']
        elif (PRD == 'ATL07'):
            VARIABLE_PATH = ['sea_ice_segments','height_segment_id']
        elif (PRD == 'ATL08'):
            VARIABLE_PATH = ['land_segments','segment_id_beg']
        elif (PRD == 'ATL10'):
            VARIABLE_PATH = ['freeboard_beam_segments','delta_time']
        elif (PRD == 'ATL12'):
            VARIABLE_PATH = ['ssh_segments','delta_time']
        # create list of valid beams within the HDF5 file
        beams = []
        for gtx in [k for k in source.keys() if bool(re.match(r'gt\d[lr]',k))]:
            # check if subsetted beam contains data
            try:
                source['/'.join([gtx,*VARIABLE_PATH])]
            except KeyError:
                pass
            else:
                beams.append(gtx)

        # for each valid beam within the HDF5 file
        frames = []
        gt = dict(gt1l=10,gt1r=20,gt2l=30,gt2r=40,gt3l=50,gt3r=60)
        for gtx in sorted(beams):
            # set variable parameters to read for specific products
            if (PRD == 'ATL06'):
                # land ice height
                var = source[gtx]['land_ice_segments']
                valid, = np.nonzero(var['h_li'][:] != var['h_li'].fillvalue)
                # variables for the output dataframe
                vnames = ['segment_id','delta_time','latitude','longitude',
                    'h_li','h_li_sigma','atl06_quality_summary',
                    'fit_statistics/dh_fit_dx',
                    'fit_statistics/dh_fit_dy',
                    'fit_statistics/dh_fit_dx_sigma',
                    'fit_statistics/n_fit_photons',
                    'fit_statistics/h_expected_rms',
                    'fit_statistics/h_robust_sprd',
                    'fit_statistics/w_surface_window_final']
            elif (PRD == 'ATL07'):
                # sea ice height
                var = source[gtx]['sea_ice_segments']
                valid, = np.nonzero(var['heights/height_segment_quality'][:] == 1)
                # variables for the output ascii file
                vnames = ['height_segment_id','seg_dist_x','delta_time',
                    'latitude','longitude',
                    'heights/height_segment_height',
                    'heights/height_segment_confidence',
                    'heights/height_segment_type',
                    'heights/height_segment_ssh_flag',
                    'heights/height_segment_w_gaussian',
                    'stats/photon_rate','stats/cloud_flag_asr',
                    'geophysical/height_segment_lpe',
                    'geophysical/height_segment_mss',
                    'geophysical/height_segment_ocean',
                    'geophysical/height_segment_ib']
            elif (PRD == 'ATL08'):
                # land and vegetation height
                var = source[gtx]['land_segments']
                valid, = np.nonzero(var['terrain/h_te_best_fit'][:] !=
                    var['terrain/h_te_best_fit'].fillvalue)
                # variables for the output dataframe
                vnames = ['segment_id_beg','segment_id_end','delta_time',
                    'latitude','longitude','brightness_flag','layer_flag',
                    'msw_flag','night_flag','terrain_flg','urban_flag',
                    'segment_landcover','segment_snowcover','segment_watermask',
                    'terrain/h_te_best_fit','terrain/h_te_uncertainty',
                    'terrain/terrain_slope','terrain/n_te_photons',
                    'canopy/h_canopy','canopy/h_canopy_uncertainty',
                    'canopy/canopy_flag','canopy/n_ca_photons']
            # create a dictionary of valid output segment values
            data = {}
            # convert data to numpy array for backwards HDF5 compatibility
            for v in vnames:
                values = np.copy(var[v][:])
                data[posixpath.basename(v)] = values[valid]
            # Generate Time Column
            delta_time = (data['delta_time']*1e9).astype('timedelta64[ns]')
            data['time'] = pandas.to_datetime(self.atlas_sdp_epoch+delta_time)
            # copy filename parameters
            data['rgt'] = np.array([int(TRK)]*len(valid))
            data['cycle'] = np.array([int(CYCL)]*len(valid))
            data['gt'] = np.array([gt[gtx]]*len(valid))
            # calculate global reference point
            if PRD in ('ATL06','ATL07','ATL08'):
                data['global_ref_pt'] = 6*1387*data[VARIABLE_PATH[-1]] + \
                    6*(data['rgt']-1) + (data['gt']/10)
            # copy beam-level attributes
            attrs = ['groundtrack_id','atlas_spot_number','atlas_beam_type',
                'sc_orientation','atmosphere_profile','atlas_pce']
            for att_name in attrs:
                att_val=self.attributes_encoder(source[gtx].attrs[att_name])
                data[att_name] = [att_val]*len(valid)
            # pandas dataframe from compiled dictionary
            frames.append(pandas.DataFrame.from_dict(data))
        # return the concatenated pandas dataframe
        return pandas.concat(frames)

    # PURPOSE: encoder for copying the file attributes
    def attributes_encoder(self, attr):
        """
        Custom encoder for copying file attributes in Python 3

        Parameters
        ----------
        attr: obj
            attribute to be converted for output
        """
        if isinstance(attr, (bytes, bytearray)):
            return attr.decode('utf-8')
        if isinstance(attr, (np.int_, np.intc, np.intp, np.int8, np.int16, np.int32,
            np.int64, np.uint8, np.uint16, np.uint32, np.uint64)):
            return int(attr)
        elif isinstance(attr, (np.float_, np.float16, np.float32, np.float64)):
            return float(attr)
        elif isinstance(attr, (np.ndarray)):
            if not isinstance(attr[0], (object)):
                return attr.tolist()
        elif isinstance(attr, (np.bool_)):
            return bool(attr)
        elif isinstance(attr, (np.void)):
            return None
        else:
            return attr
