# stdlib imports
from collections import OrderedDict

# third party imports
from impactutils.io.table import read_excel

UNITS = {'PGA': r'%g',
         'PGV': 'cm/s',
         'SA(0.3)': r'%g',
         'SA(1.0)': r'%g',
         'SA(3.0)': r'%g',
         }


def get_station_feature(trow):
    index = trow.index.to_list()
    row = {}
    components = {}
    for idx_tpl in index:
        if isinstance(idx_tpl, str):
            row[idx_tpl] = trow[idx_tpl]
        else:
            if 'UNNAMED' in idx_tpl[1]:
                row[idx_tpl[0]] = trow[idx_tpl[0]][idx_tpl[1]]
            else:
                component = idx_tpl[0]
                imt = idx_tpl[1]
                if component not in components:
                    components[component] = {imt: trow[component][imt]}
                else:
                    components[component][imt] = trow[component][imt]

    scode = '%s.%s' % (row['NETID'], row['STATION'])
    station_feature = OrderedDict()
    station_properties = OrderedDict()
    station_feature['type'] = 'Feature'
    station_feature['id'] = scode
    station_properties['name'] = ''
    if 'LOC' in row:
        station_properties['name'] = row['LOC']

    station_properties['code'] = row['STATION']
    if 'INTENSITY' in row:
        station_properties['network'] = 'INTENSITY'
        # limit MMI to between 1-10
        if row['INTENSITY'] > 10.0:
            row['INTENSITY'] = 10.0
        if row['INTENSITY'] < 1.0:
            row['INTENSITY'] = 1.0
        station_properties['intensity'] = row['INTENSITY']
        station_properties['intensity_flag'] = row['FLAG']
        station_properties['intensity_stddev'] = 0.0
        station_properties['station_type'] = 'macroseismic'
    else:
        station_properties['network'] = row['NETID']
        # station_properties['station_type'] = 'seismic'
    if 'DISTANCE' in row:
        station_properties['distance'] = row['DISTANCE']

    # station_properties['source'] = stream[0].stats.standard['source']
    station_properties['source'] = row['SOURCE']

    coordinates = (row['LON'], row['LAT'])
    station_feature['geometry'] = {'type': 'Point',
                                   'coordinates': coordinates}

    if not len(components):
        station_properties['channels'] = []
    else:
        # get the channels from the row index
        channels = []
        for component, imtdict in components.items():
            compdict = {}
            compdict['name'] = component
            amps = []
            for imtname, imtvalue in imtdict.items():
                ampdict = {}
                ampdict['name'] = imtname.lower()
                ampdict['value'] = imtvalue
                ampdict['units'] = UNITS[imtname]
                ampdict['flag'] = 0
                ampdict['ln_sigma'] = 0
                amps.append(ampdict)
            compdict['amplitudes'] = amps
            channels.append(compdict)
        station_properties['channels'] = channels

    station_feature['properties'] = station_properties
    return station_feature


def get_stations_dict(excelfile):
    dataframe, reference = read_excel(excelfile)

    features = []
    for idx, row in dataframe.iterrows():
        feature = get_station_feature(row)
        features.append(feature)

    jdict = {'type': 'FeatureCollection', 'features': features}
    return (jdict, len(dataframe))
