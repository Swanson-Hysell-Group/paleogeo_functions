# standard modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.cm import get_cmap
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

# pmagpy
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import cartopy.crs as ccrs
import cartopy
import xml.etree.ElementTree as ET
import pygplates as pgp

def create_vgp_FeatureCollection(compilation):
    """
    Loop through poles and produce a pygplates FeatureCollection of multiple VGPs.
    
    Modified from code by Michael G. Tetley.
    
    Parameters
    ----------
    compilation : dataframe
        pole compilation
        
    Returns
    -------
    vpgFeatureCollection : FeatureCollection
        pygplates FeatureCollection of VGPs in compilation
    """
    poleLat = []
    poleLon = []
    poleName = []
    poleSiteLat = []
    poleSiteLon = []
    poleNominalAge = []
    poleA95 = []
    poleAgeLowerLimit = []
    poleAgeUpperLimit = []
    plateID = []

    count = 0

    for i in range(len(compilation)):

        if np.isfinite(compilation['slat'][i]) and np.isfinite(compilation['slon'][i]) and \
           np.isfinite(compilation['age lower'][i]) and np.isfinite(compilation['age upper'][i]):

            poleLat.append(compilation['plat'][i])
            poleLon.append(compilation['plon'][i])

            poleName.append(compilation['name'][i] + ' (' + compilation['grade'][i] + ')')
            poleSiteLat.append(compilation['slat'][i])
            poleSiteLon.append(compilation['slon'][i])
            poleNominalAge.append(compilation['age'][i])
            poleA95.append(compilation['a95'][i])

            poleAgeLowerLimit.append(compilation['age lower'][i])
            poleAgeUpperLimit.append(compilation['age upper'][i])

            plateID.append(compilation['plateID'][i])

            count = count + 1

        # Print if any of the isfinite tests fail
        else:

            print('Bad data for : {}'.format(compilation['name'][i]))


    # Create new GPlates Feature Collection
    vpgFeatureCollection = pgp.FeatureCollection()

    # Create new GPlates feature 'VirtualGeomagneticPole'.
    # Pole lat, pole lon, pole name, and reconstruction plate ID added within PointOnSphere method.
    # Inc, Dec, A95, Age and Sample site lat/lon values to added within 'other_properties' method.

    for j in range(count):

        vgpFeature = pgp.Feature.create_reconstructable_feature(
                     pgp.FeatureType.create_gpml('VirtualGeomagneticPole'),
                     pgp.PointOnSphere([np.float(poleLat[j]), np.float(poleLon[j])]),
                     name = poleName[j],
                     reconstruction_plate_id = int(plateID[j]),
                     other_properties = [(pgp.PropertyName.create_gpml('poleA95'),
                                          pgp.XsDouble(np.float64(poleA95[j]))),
                                         (pgp.PropertyName.create_gpml('averageAge'),
                                          pgp.XsDouble(np.float64(poleNominalAge[j]))),
                                         (pgp.PropertyName.create_gpml('averageSampleSitePosition'),
                                          pgp.GmlPoint(pgp.PointOnSphere([np.float(poleSiteLat[j]), 
                                                                          np.float(poleSiteLon[j])])))])

        # Add newly created feature to existing Feature Collection
        vpgFeatureCollection.add(vgpFeature)

    return vpgFeatureCollection


def create_slim_vgp_FeatureCollection(compilation):
    """
    Loop through poles and produce a pygplates FeatureCollection of multiple VGPs.
    
    Modified from code by Michael G. Tetley.
    
    Parameters
    ----------
    compilation : dataframe
        pole compilation
        
    Returns
    -------
    vpgFeatureCollection : FeatureCollection
        pygplates FeatureCollection of VGPs in compilation
    """
    poleLat = []
    poleLon = []
    poleName = []
    poleNominalAge = []
    poleA95 = []
    plateID = []

    count = 0

    for i in range(len(compilation)):

        poleLat.append(compilation['plat'][i])
        poleLon.append(compilation['plon'][i])

        poleName.append(compilation['name'][i] + ' (' + str(compilation['grade'][i]) + ')')
        poleNominalAge.append(compilation['age'][i])
        poleA95.append(compilation['a95'][i])

        plateID.append(compilation['plateID'][i])

        count = count + 1

    # Create new GPlates Feature Collection
    vpgFeatureCollection = pgp.FeatureCollection()

    # Create new GPlates feature 'VirtualGeomagneticPole'.
    # Pole lat, pole lon, pole name, and reconstruction plate ID added within PointOnSphere method.
    # Inc, Dec, A95, Age and Sample site lat/lon values to added within 'other_properties' method.

    for j in range(count):

        vgpFeature = pgp.Feature.create_reconstructable_feature(
                     pgp.FeatureType.create_gpml('VirtualGeomagneticPole'),
                     pgp.PointOnSphere([np.float(poleLat[j]), np.float(poleLon[j])]),
                     name = poleName[j],
                     reconstruction_plate_id = int(plateID[j]),
                     other_properties = [(pgp.PropertyName.create_gpml('poleA95'), pgp.XsDouble(poleA95[j])),
                                         (pgp.PropertyName.create_gpml('averageAge'), pgp.XsDouble(poleNominalAge[j]))])

        # Add newly created feature to existing Feature Collection
        vpgFeatureCollection.add(vgpFeature)

    return vpgFeatureCollection

    
def create_vgp_gpml(vpgFeatureCollection, filename):
    """
    Create a .gpml for a FeatureCollection of VGPs.
    
    Modified from code by Michael G. Tetley.
    
    Parameters
    ----------
    vpgFeatureCollection : FeatureCollection
        pygplates FeatureCollection of VGPs in compilation
        
    filename : string
        path and name for output .gpml
    """
    # Generate GPML output file
    gpmlOutputFile = filename

    # Check for existing output file with same name and remove if found
    if os.path.isfile(filename):
        os.remove(filename)

    # Check to make sure vgpFeatureCollection (feature collection) is not empty before writing to file
    if len(vpgFeatureCollection) != 0:
        outputFeatureCollection = pgp.FeatureCollectionFileFormatRegistry()
        outputFeatureCollection.write(vpgFeatureCollection, filename)

    # Check if new file was created and confirm export
    if os.path.isfile(filename):
        print('Palaeomagnetic pole data successfully exported in GPML format.')
        
        
def get_craton_XYs(gpml, plateIDs):
    """
    Get XY coordinates of a plate polygon from a .gpml.
    
    Parameters
    ----------
    gpml : string
        Path to .gpml file.
        
    plateIDs : list
        Of plateIDs.
    """
    # namespace dictionary
    ns = {'gpml':'http://www.gplates.org/gplates',
          'gml':'http://www.opengis.net/gml'}
    
    # initial parse
    tree = ET.parse(gpml)
    root = tree.getroot()
    
    # storage
    Xs = []
    Ys = []
    
    # iterate through featureMembers
    for featureMember in root.findall('gml:featureMember',ns):
        
        # get child
        for child in featureMember:
            slice_ind = child.tag.find('}')
            child_root = 'gpml:' + child.tag[slice_ind+1:]
        
        # check plateID
        plateID_path = child_root + '/gpml:reconstructionPlateId/gpml:ConstantValue/gpml:value'
        feature_plateID = int(featureMember.find(plateID_path,ns).text)
        if feature_plateID in plateIDs:
            
            if featureMember.find(child_root + '/gpml:outlineOf', ns)!=None:
                polygon_root = child_root + '/gpml:outlineOf'
            elif featureMember.find(child_root + '/gpml:boundary', ns)!=None:
                polygon_root = child_root + '/gpml:boundary'
            elif featureMember.find(child_root + '/gpml:unclassifiedGeometry', ns)!=None:
                polygon_root = child_root + '/gpml:unclassifiedGeometry'
            elif featureMember.find(child_root + '/gpml:centerLineOf', ns)!=None:
                polygon_root = child_root + '/gpml:centerLineOf'
            else:
                raise Exception('polygon_root undefined.')
            
            # get coordinates
            posList_path = polygon_root + '/gpml:ConstantValue/gpml:value/gml:Polygon/gml:exterior/gml:LinearRing/gml:posList'
            for feature_posList in featureMember.findall(posList_path,ns):
                np_posList = np.fromstring(feature_posList.text, dtype=float, sep=' ')
            
                # split into lat and lon
                lat_inds = np.arange(0, len(np_posList), 2, dtype=int)
                lon_inds = np.arange(1, len(np_posList), 2, dtype=int)

                feature_lat = np_posList[lat_inds]
                feature_lon = np_posList[lon_inds]
            
                Xs.append(feature_lon)
                Ys.append(feature_lat)
            
    return Xs, Ys


def get_single_craton_XYs(gpml):
    """
    Get XY coordinates of a plate polygon from a .gpml.
    
    Parameters
    ----------
    gpml : string
        Path to .gpml file.
    """
    # namespace dictionary
    ns = {'gpml':'http://www.gplates.org/gplates',
          'gml':'http://www.opengis.net/gml'}
    
    # initial parse
    tree = ET.parse(gpml)
    root = tree.getroot()
    
    # storage
    Xs = []
    Ys = []
    
    # iterate through featureMembers
    featureMember = root.find('gml:featureMember',ns)
        
    # get child
    for child in featureMember:
        slice_ind = child.tag.find('}')
        child_root = 'gpml:' + child.tag[slice_ind+1:]

    if featureMember.find(child_root + '/gpml:outlineOf', ns)!=None:
        polygon_root = child_root + '/gpml:outlineOf'
    elif featureMember.find(child_root + '/gpml:boundary', ns)!=None:
        polygon_root = child_root + '/gpml:boundary'
    elif featureMember.find(child_root + '/gpml:unclassifiedGeometry', ns)!=None:
        polygon_root = child_root + '/gpml:unclassifiedGeometry'
    elif featureMember.find(child_root + '/gpml:centerLineOf', ns)!=None:
        polygon_root = child_root + '/gpml:centerLineOf'
    else:
        raise Exception('polygon_root undefined.')

    # get coordinates
    posList_path = polygon_root + '/gpml:ConstantValue/gpml:value/gml:Polygon/gml:exterior/gml:LinearRing/gml:posList'
    for feature_posList in featureMember.findall(posList_path,ns):
        np_posList = np.fromstring(feature_posList.text, dtype=float, sep=' ')

        # split into lat and lon
        lat_inds = np.arange(0, len(np_posList), 2, dtype=int)
        lon_inds = np.arange(1, len(np_posList), 2, dtype=int)

        feature_lat = np_posList[lat_inds]
        feature_lon = np_posList[lon_inds]

        Xs = feature_lon
        Ys = feature_lat
            
    return Xs, Ys



def pt_rot_vec(EP, Lats, Lons):
    """
    Vectorized Euler pole rotation (Cox & Hart 1986, Box 7-3).
    
    Parameters
    ----------
    EP : [lat, lon, angle]
        Euler pole latitude, longitude (degrees), rotation angle (degrees, CCW).
    Lats : array-like
        Latitudes of points (deg).
    Lons : array-like
        Longitudes of points (deg).
    
    Returns
    -------
    RLats : ndarray
        Rotated latitudes (deg).
    RLons : ndarray
        Rotated longitudes (deg).
    """

    # ensure numpy arrays
    Lats = np.asarray(Lats, dtype=float)
    Lons = np.asarray(Lons, dtype=float)

    # build rotation matrix once
    E = pmag.dir2cart([EP[1], EP[0], 1.])  # pole direction (lon, lat, radius=1)
    omega = np.radians(EP[2])
    c, s = np.cos(omega), np.sin(omega)
    Ex, Ey, Ez = E

    R = np.array([
        [Ex*Ex*(1-c)+c,     Ex*Ey*(1-c)-Ez*s, Ex*Ez*(1-c)+Ey*s],
        [Ey*Ex*(1-c)+Ez*s, Ey*Ey*(1-c)+c,     Ey*Ez*(1-c)-Ex*s],
        [Ez*Ex*(1-c)-Ey*s, Ez*Ey*(1-c)+Ex*s, Ez*Ez*(1-c)+c]
    ])

    # mask for delimiters (Lats > 90)
    mask = Lats <= 90

    # convert input lat/lon to Cartesian
    x, y, z = pmag.dir2cart(np.array([Lons[mask], Lats[mask], np.ones_like(Lats[mask])]).T).T

    # apply rotation (matrix multiply on all points)
    A = np.vstack([x, y, z])              # shape (3, N)
    Ap = R @ A                            # shape (3, N)

    # convert back to lon/lat
    lon_rot, lat_rot, _ = pmag.cart2dir(np.array([Ap[0], Ap[1], Ap[2]]).T).T

    # prepare outputs
    RLats = np.full_like(Lats, np.nan)
    RLons = np.full_like(Lons, np.nan)
    RLats[mask] = lat_rot
    RLons[mask] = lon_rot
    RLats[~mask] = Lats[~mask]  # preserve delimiters
    RLons[~mask] = Lons[~mask]

    return RLats, RLons

def craton_plot(ax, plateIDs, Eulers, edgecolor, facecolor, alpha, linewidth, gpml = '../GPlates/Cratons/shapes_cratons.gpml', reverse_draw=False,draw_face=True,draw_edge=True):
    """
    Plot cratons with rotation.
    
    Parameters
    ----------
    ax : map axis
        On which to plot.
    
    plateIDs : list
        Of plateIDs.
        
    Eulers : list of lists
        Of Euler rotation parameters - if more than one given,
        the rotations will be additive.
    """
    # get cratons from .gpml
    
    Xs, Ys = get_craton_XYs(gpml, plateIDs)
    
    # draw in reverse
    if reverse_draw:
        Xs = np.flip(Xs)
        Ys = np.flip(Ys)
    
    # rotate cratons
    rotated_Xs = []
    rotated_Ys = []
    for i in range(len(Xs)):
        x = np.array(Xs[i])
        y = np.array(Ys[i])
        for euler in Eulers:
            y, x = pmag.pt_rot(euler, y, x)   # ideally accepts vectors
        rotated_Xs.append(x)
        rotated_Ys.append(y)
        
    # add cratons

    patches_list = []
    for x, y in zip(rotated_Xs, rotated_Ys):
        XY = np.stack([x[::-1], y[::-1]], axis=1)
        patches_list.append(patches.Polygon(XY, closed=True))

    # one collection for edges
    if draw_edge:
        edge_collection = PatchCollection(patches_list, facecolor='none',
                                        edgecolor=edgecolor, alpha=alpha,
                                        linewidth=linewidth, transform=ccrs.Geodetic())
        ax.add_collection(edge_collection)

    # one collection for faces
    if draw_face:
        face_collection = PatchCollection(patches_list, facecolor=facecolor,
                                        edgecolor='none', alpha=alpha,
                                        transform=ccrs.Geodetic())
        ax.add_collection(face_collection)

        
        
def single_craton_plot(ax, gpml, Eulers, edgecolor, facecolor, alpha, linewidth):
    """
    Plot cratons with rotation.
    
    Parameters
    ----------
    ax : map axis
        On which to plot.
    
    gpml : string
        Path to .gpml file.
        
    Eulers : list of lists
        Of Euler rotation parameters - if more than one given,
        the rotations will be additive.
    """
    # get cratons from .gpml
    Xs, Ys = get_single_craton_XYs(gpml)
    
    # rotate craton
    rotated_Xs = np.array([])
    rotated_Ys = np.array([])
    for i in range(len(Xs)):
        this_X = [Xs[i]]
        this_Y = [Ys[i]]
        for j in range(len(Eulers)):
            this_Y, this_X = pmag.pt_rot(Eulers[j], this_Y, this_X)
        rotated_Xs = np.append(rotated_Xs, this_X)
        rotated_Ys = np.append(rotated_Ys, this_Y)
        
    # add craton
    #XY = np.stack([rotated_Xs[::-1],rotated_Ys[::-1]],axis=1)
    XY = np.stack([rotated_Xs,rotated_Ys],axis=1)
    print(XY.shape)
    poly_edge = patches.Polygon(XY,
                                edgecolor=edgecolor,facecolor='none',alpha=alpha,
                                transform=ccrs.Geodetic(), linewidth=linewidth)
    poly_face = patches.Polygon(XY,
                                edgecolor='none',facecolor=facecolor,alpha=alpha,
                                transform=ccrs.Geodetic())
    ax.add_patch(poly_face)
    ax.add_patch(poly_edge)
    
    
def equi_filled(map_axis, centerlon, centerlat, radius, color, alpha=1.0, edge_alpha=1.0):
    """
    Modified from the ipmag function equi().
    """
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = ipmag.shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
    
    X = X[::-1]
    Y = Y[::-1]
    
    XY = np.stack([X,Y],axis=1)
    
    circle_edge = patches.Polygon(XY,
                                  edgecolor=color,facecolor='none',alpha=edge_alpha,
                                  transform=ccrs.Geodetic())
    circle_face = patches.Polygon(XY,
                                  edgecolor='none',facecolor=color,alpha=alpha,
                                  transform=ccrs.Geodetic())
    
    map_axis.add_patch(circle_face)
    map_axis.add_patch(circle_edge)
    
    
def rotated_pole_plot(ax, plon, plat, a95, Eulers, marker, s, marker_color, a95_color, a95_alpha, a95_edge_alpha=1.0):
    """
    Plot paleomagnetic pole with rotation.
    """
    # rotate pole
    rotated_plat = plat
    rotated_plon = plon
    for i in range(len(Eulers)):
        rotated_plat, rotated_plon = pmag.pt_rot(Eulers[i], [rotated_plat], [rotated_plon])
        rotated_plat = rotated_plat[0]
        rotated_plon = rotated_plon[0]
    
    # degrees to km conversion
    a95_km = a95 * 111.32
    
    # pole
    ax.scatter(rotated_plon, rotated_plat, marker=marker,
               color=marker_color, edgecolors='k', s=s,
               label='__nolegend__', zorder=101, transform=ccrs.PlateCarree())
    
    # a95
    equi_filled(ax, rotated_plon, rotated_plat, a95_km, a95_color, alpha=a95_alpha, edge_alpha=a95_edge_alpha)
    

def rotated_point_plot(ax, plon, plat, Eulers, s, marker_color):
    """
    Plot point with rotation.
    """
    # rotate pole
    rotated_plat = plat
    rotated_plon = plon
    for i in range(len(Eulers)):
        rotated_plat, rotated_plon = pmag.pt_rot(Eulers[i], rotated_plat, rotated_plon)
        rotated_plat = rotated_plat
        rotated_plon = rotated_plon
    
    # pole
    ax.scatter(rotated_plon, rotated_plat,
               color=marker_color, s=s, edgecolor='none',
               label='__nolegend__', zorder=101, transform=ccrs.PlateCarree())
    
def APWP_plot(ax,
              plon_rot, plat_rot, a95_rot, age_rot,
              plon_fix, plat_fix, a95_fix, age_fix,
              Euler,
              label_rot, color_rot, marker_rot, s_rot,
              label_fix, color_fix, marker_fix, s_fix,
              age_lim=None, cmap='viridis'):
    """
    Plot apparent polar wander paths for two connected cratons.
    """
    # get vmin and vmax, and slice data if necessary
    if age_lim==None:
        vmin = np.min([age_rot.min(), age_fix.min()])
        vmax = np.max([age_rot.max(), age_fix.max()])
    else:
        vmin = age_lim[0]
        vmax = age_lim[1]
        
        rot_mask = (age_rot>=age_lim[0])&(age_rot<=age_lim[1])
        plon_rot = plon_rot[rot_mask]
        plat_rot = plat_rot[rot_mask]
        a95_rot = a95_rot[rot_mask]
        age_rot = age_rot[rot_mask]
        
        fix_mask = (age_fix>=age_lim[0])&(age_fix<=age_lim[1])
        plon_fix = plon_fix[fix_mask]
        plat_fix = plat_fix[fix_mask]
        a95_fix = a95_fix[fix_mask]
        age_fix = age_fix[fix_mask]
    
    # rotate poles
    rotated_plat = np.array([])
    rotated_plon = np.array([])
    
    for i in range(len(plon_rot)):
        this_plat, this_plon = pmag.pt_rot(Euler, [plat_rot[i]], [plon_rot[i]])
        rotated_plat = np.append(rotated_plat, this_plat[0])
        rotated_plon = np.append(rotated_plon, this_plon[0])
    
    # degrees to km conversion
    a95_rot_km = a95_rot * 111.32
    a95_fix_km = a95_fix * 111.32
    
    # colormap to age
    color_mapping = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    colors_rot = color_mapping.to_rgba(age_rot).tolist()
    colors_fix = color_mapping.to_rgba(age_fix).tolist()
    
    # pole
    ax.scatter(plon_fix, plat_fix, marker=marker_fix,
               color=colors_fix, edgecolors='k', s=s_fix,
               label='__nolegend__', zorder=101, transform=ccrs.PlateCarree())
    ax.scatter(rotated_plon, rotated_plat, marker=marker_rot,
               color=colors_rot, edgecolors='k', s=s_rot,
               label='__nolegend__', zorder=101, transform=ccrs.PlateCarree())
    
    # a95
    for i in range(len(plon_fix)):
        equi_filled(ax, plon_fix[i], plat_fix[i], a95_fix_km[i], color_fix, alpha=0.3)
    for i in range(len(rotated_plon)):
        equi_filled(ax, rotated_plon[i], rotated_plat[i], a95_rot_km[i], color_rot, alpha=0.3)
        
    # create fake legend
    ax.scatter([], [], marker=marker_fix,
               color=color_fix, edgecolors='k', s=s_fix,
               label=label_fix, transform=ccrs.PlateCarree())
    ax.scatter([], [], marker=marker_rot,
               color=color_rot, edgecolors='k', s=s_rot,
               label=label_rot, transform=ccrs.PlateCarree())
    
    # colorbar
    color_mapping._A = []
    plt.colorbar(color_mapping, orientation='horizontal', shrink=0.8,
                 pad=0.05, label='age [Ma]')
    
    # prettify
    ax.legend(loc=4, fontsize=12)
    ax.set_title('{} rotated to {} - Euler : {}'.format(label_rot,label_fix,str(Euler)))
    
    plt.show()
    
    
def motion_path_calc(point_lat, point_lon, min_age, max_age, rotation_file_path, moving_plate_id, relative_plate_id):
    """
    Calculate the motion path of a point given a rotation model.
    
    Parameters
    ----------
    point_lat : float
        Latitude of point to follow.
        
    point_lon : float
        Longitude of point to follow.
        
    min_age : float
        Age at which to stop the reconstruction (Ma).
        
    max_age : float
        Age at which to begin the reconstruction (Ma).
        
    rotation_file_path : string
        Path to rotation file.
        
    moving_plate_id : int
        Plate ID to which the point belongs.
        
    relative_plate_id : int
        Plate ID of the relative plate.
        
    Returns
    -------
    times_output : array
        Age array.
        
    trail : array
        With (lat, lon) of the point.
    """
    min_age = np.float64(min_age)
    max_age = np.float64(max_age)
    
    rotation_model = pgp.RotationModel(rotation_file_path)
    
    SeedPoint = (point_lat,point_lon)
    MovingPlate = moving_plate_id
    RelativePlate = relative_plate_id
    times = np.arange(0,max_age,2)
    times_output = np.arange(min_age,max_age,2)
    
    # Create a motion path feature
    digitisation_time = 0
    seed_points_at_digitisation_time = pgp.MultiPointOnSphere([SeedPoint]) 
    motion_path_feature = pgp.Feature.create_motion_path(
            seed_points_at_digitisation_time,
            times,
            valid_time=(max_age,0),
            relative_plate=RelativePlate,
            reconstruction_plate_id = MovingPlate)
    
    # Create the shape of the motion path
    reconstruction_time = 0
    reconstructed_motion_paths = []
    pgp.reconstruct(
            motion_path_feature, rotation_model, reconstructed_motion_paths, reconstruction_time,
            reconstruct_type=pgp.ReconstructType.motion_path)

    # get the reconstructed coordinates into numpy arrays
    for reconstructed_motion_path in reconstructed_motion_paths:
        trail = reconstructed_motion_path.get_motion_path().to_lat_lon_array()
        
    times_output = np.flipud(times_output)
    trail = trail[0:len(times_output)]
        
    return times_output, trail


def centroid_analysis(centroid_lat, centroid_lon, df, 
                      TCedit_plateID, M21_plateID,
                      TC_edit_rot_path = 'models/extended_TC17/TC2017-SHM2017-D2018-extended.rot',
                      M21_rot_path = 'models/M21/1000_0_rotfile_Merdith_et_al.rot',
                      age_max_TC=1150,age_min_TC=250,
                      age_max_M21=1000,age_min_M21=250,):
    """
    Perform a series of calculations for a centroid point on a given craton.
    
    Parameters
    ----------
    centroid_lat : float
        Centroid latitude.
    centroid_lon : float
        Centroid longitude.
    df : dataframe
        With pole data. 
    TCedit_plateID : int
        plateID for the Torsvik and Cocks (2012) - edited model. 
    M21_plateID : int
        plateID for the Merdith et al. (2020) model.
    age_max : int
        maximum age of analysis and plot
    age_min : int
        minimum age of analysis and plot
        
    Returns
    -------
    df : dataframe
        The pole data, but with some additional columns.
        
    TC_edited_centroid_path : dataframe
        Centroid path implied by the Torsvik and Cocks (2012) - edited model.
        
    M21_centroid_path : dataframe
        Centroid path implied by the Merdith et al. (2020) model.
    """
    
    # calculate implied centroid paleolatitude
    df['centroid_lat'] = pd.Series()

    for i in range(len(df)):
        lat = ipmag.lat_from_pole(centroid_lon, centroid_lat,
                                  df['plon'][i], df['plat'][i])
        df['centroid_lat'][i] = lat
        
    # calculate centroid motion based on .rot files
    TC_edited_times, TC_edited_trail = motion_path_calc(centroid_lat, centroid_lon,
                                                  age_min_TC, age_max_TC, TC_edit_rot_path, TCedit_plateID, 1)
    TC_edited_centroid_path = pd.DataFrame({'age':TC_edited_times, 'lat':TC_edited_trail[:,0], 'lon':TC_edited_trail[:,1]})

    M21_times, M21_trail = motion_path_calc(centroid_lat, centroid_lon,
                                                  age_min_M21, age_max_M21, M21_rot_path, M21_plateID, 0)
    M21_centroid_path = pd.DataFrame({'age':M21_times, 'lat':M21_trail[:,0], 'lon':M21_trail[:,1]})
    
    # create +/- age columns
    df['age upper diff'] = df['age upper'] - df['age']
    df['age lower diff'] = df['age'] - df['age lower']
    
    # plot
    df_A = df[df['grade']=='A']
    df_B = df[df['grade']=='B']
    df_C = df[df['grade']=='C+']
    df_NR = df[df['grade']=='NR']

    cmap = get_cmap('Blues')
    A_pole_color = 'red'
    B_pole_color = cmap(0.8)
    C_pole_color = cmap(0.2)
    NR_pole_color = cmap(0.5)

    fig, ax = plt.subplots(figsize=(10,5))

    ax.errorbar(df_A['age'], df_A['centroid_lat'],
                yerr=df_A['a95'], xerr=[df_A['age lower diff'], df_A['age upper diff']],
                fmt='o', color=A_pole_color, label='A poles')

    ax.errorbar(df_B['age'], df_B['centroid_lat'],
                yerr=df_B['a95'], xerr=[df_B['age lower diff'], df_B['age upper diff']],
                fmt='o', color=B_pole_color, label='B poles')
    
    ax.errorbar(df_C['age'], df_C['centroid_lat'],
                yerr=df_C['a95'], xerr=[df_C['age lower diff'], df_C['age upper diff']],
                fmt='o', color=C_pole_color, label='C poles')

    ax.errorbar(df_NR['age'], df_NR['centroid_lat'],
                yerr=df_NR['a95'], xerr=[df_NR['age lower diff'], df_NR['age upper diff']],
                fmt='o', color=NR_pole_color, label='NR poles')

    ax.plot(M21_centroid_path['age'], M21_centroid_path['lat'],
            c='C7', lw=4, zorder=-90, alpha=0.5, label='Merdith et al. (2021)')

    ax.plot(TC_edited_centroid_path['age'], TC_edited_centroid_path['lat'],
            c='C7', lw=4, zorder=-99, alpha=0.9, label='prelim model + Torsvik and Cocks (2017)')

    ax.set_ylim(-90,90)
    ax.set_xlim(age_max_TC, age_min_TC)
    ax.set_xlabel('age [Ma]')
    ax.set_ylabel('paleolatitude [$^{\circ}$]')
    ax.legend()

    plt.show(fig)
    
    # returns
    return df, TC_edited_centroid_path, M21_centroid_path

def find_plate_ids_in_region(feature_path, id_prefix, time):
    '''
    Find all plate IDs in a given .gpml feature collection that:
      - Start with a specific prefix
      - Have a begin time older than the time of interest

    Parameters:
    feature_path : str
        Path to the .gpml feature collection file.
    id_prefix : str or int
        Prefix of the plate IDs to search for.
    time : float
        Time of interest (in Ma).

    Returns:
    list
        Sorted list of unique plate IDs that match the criteria.
    '''
    features = pgp.FeatureCollection(feature_path)
    plate_ids = []
    for feature in features:
        pid = feature.get_reconstruction_plate_id()
        if str(pid).startswith(str(id_prefix)):
            begin_time, end_time = feature.get_valid_time()
            if begin_time >= time:   # begin time older than time of interest
                plate_ids.append(pid)
    return sorted(set(plate_ids))   # unique & sorted

def plot_region(ax, feature_path, rotation_model, region_id, fixed_plate, time, color='lightgrey', cratons_alpha=0.65, lw=0.5, edgecolor='k'):
    '''
    Plot a region on the map given a feature collection and rotation model.
    The region is defined by the region_id prefix, and all plates with IDs starting
    with that prefix and valid at the given time are plotted.

    Parameters:
    ax : matplotlib.axes.Axes
        The axes on which to plot.
    feature_path : str
        Path to the .gpml feature collection file.
    rotation_model : pgp.RotationModel
        The rotation model to use for plate rotations.
    region_id : str or int
        Prefix of the plate IDs defining the region.
    fixed_plate : int
        The fixed plate ID for the rotation model.
    time : float
        Time of interest (in Ma).
    color : str, optional
        Fill color for the region (default is 'lightgrey').
    cratons_alpha : float, optional
        Transparency for the region fill (default is 0.65).
    lw : float, optional
        Line width for the region outline (default is 0.5).
    edgecolor : str, optional
        Edge color for the region outline (default is 'k').
    '''
    plate_ids = find_plate_ids_in_region(feature_path, region_id, time)
    for pid in plate_ids:
        rotation = rotation_model.get_rotation(time, pid, 0, fixed_plate, 1).get_lat_lon_euler_pole_and_angle_degrees()
        craton_plot(ax, [pid], [rotation],
                    edgecolor, color, cratons_alpha, lw, gpml=feature_path, reverse_draw=False)
        
def plot_craton(ax, feature_path, rotation_model, plate_id, fixed_plate, time, color='lightgrey', cratons_alpha=0.65, lw=0.5, edgecolor='k'):
    '''
    Plot a single craton on the map given a feature collection and rotation model.

    Parameters:
    ax : matplotlib.axes.Axes
        The axes on which to plot.
    feature_path : str
        Path to the .gpml feature collection file.
    rotation_model : pgp.RotationModel
        The rotation model to use for plate rotations.
    plate_id : int
        The plate ID of the craton to plot.
    fixed_plate : int
        The fixed plate ID for the rotation model.
    time : float
        Time of interest (in Ma).
    color : str, optional
        Fill color for the craton (default is 'lightgrey').
    cratons_alpha : float, optional
        Transparency for the craton fill (default is 0.65).
    lw : float, optional
        Line width for the craton outline (default is 0.5).
    edgecolor : str, optional
        Edge color for the craton outline (default is 'k').
    '''
    rotation = rotation_model.get_rotation(time, plate_id, 0, fixed_plate, 1).get_lat_lon_euler_pole_and_angle_degrees()
    craton_plot(ax, [plate_id], [rotation],
                edgecolor, color, cratons_alpha, lw, gpml=feature_path, reverse_draw=False)
    

def plot_rotated_craton(ax, feature_path, rotation_model, plate_id, fixed_plate, time, euler, angle, color='lightgrey', cratons_alpha=0.65, lw=0.5, edgecolor='k'):
    '''
    Plot a single craton on the map given a feature collection and rotation model, with an additional arbitrary rotation applied.
    This function is a helper function for testing rotation of a craton by a specified Euler pole and angle.

    Parameters:
    ax : matplotlib.axes.Axes
        The axes on which to plot.
    feature_path : str
        Path to the .gpml feature collection file.
    rotation_model : pgp.RotationModel
        The rotation model to use for plate rotations.
    plate_id : int
        The plate ID of the craton to plot.
    fixed_plate : int
        The fixed plate ID for the rotation model.
    time : float
        Time of interest (in Ma).
    euler : tuple
        Euler pole coordinates (latitude, longitude) for the additional rotation.
    angle : float
        Rotation angle (in degrees) for the additional rotation.
    color : str, optional
        Fill color for the craton (default is 'lightgrey').
    cratons_alpha : float, optional
        Transparency for the craton fill (default is 0.65).
    lw : float, optional
        Line width for the craton outline (default is 0.5).
    edgecolor : str, optional
        Edge color for the craton outline (default is 'k').
    '''
    rotation = rotation_model.get_rotation(time, plate_id, 0, fixed_plate, 1)
    adjust_euler = pgp.FiniteRotation(euler, np.deg2rad(angle))
    test_finite_rotation = adjust_euler * rotation
    craton_plot(ax, [plate_id], [test_finite_rotation.get_lat_lon_euler_pole_and_angle_degrees()],
                edgecolor, color, cratons_alpha, lw, gpml=feature_path, reverse_draw=False)
    