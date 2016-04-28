#!/usr/bin/env python

#stdlib imports
import os.path
import sys
import io
import tempfile
import textwrap

#third party
import numpy as np

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..','..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

from shakemap.mapping.gmtcolormap import GMTColorMap


def test():
    cptfiletext1 = '''
    0       255     255     255     1       255     255     255
    1       255     255     255     2       191     204     255
    2       191     204     255     3       160     230     255
    3       160     230     255     4       128     255     255
    4       128     255     255     5       122     255     147
    5       122     255     147     6       255     255     0
    6       255     255     0       7       255     200     0
    7       255     200     0       8       255     145     0
    8       255     145     0       9       255     0       0
    9       255     0       0       10      200     0       0
    #10      200     0       0       13      128     0       0
    '''
    cptfiletext2 = '''
    #COLOR_MODEL = RGB
    #
    -100    195     255     193     0       110     135     80
    0       110     135     80      50      120     160     90
    50      120     160     90      350     230     220     110
    350     230     220     110     1000    210     170     80
    1000    210     170     80      1800    195     140     100
    1800    195     140     100     2300    100     80      70
    2300    100     80      70      2600    60      60      60
    2600    60      60      60      4000    255     255     255
    4000    255     255     255     9000    255     255     255
    9000    255     255     255     9100    255     128     0
    9100    255     128     0       9200    255     0       0
    B       0       0       0
    F       255     255     255
    N       128     128     128
    '''
    cptfiles = [cptfiletext1,cptfiletext2]
    vmins = [0,-100]
    vmaxs = [10,9200]
    vmids = [5,2600]
    cmins = [(1.0,1.0,1.0),([ 0.76470588, 1.0, 0.75686275])]
    cmaxs = [(0.78431373, 0.0, 0.0),(1.0,0.0,0.0)]
    cmids = [(0.48865821, 1.0, 0.56516724),(0.23590927, 0.23560169, 0.2354479)]
    midnorms = [0.5,0.29032258]
    try:
        for i in range(0,len(cptfiles)):
            tmp,tfile = tempfile.mkstemp()
            os.close(tmp)
            f = open(tfile,'wt')
            cptfiletext = cptfiles[i]
            vmin = vmins[i]
            vmid = vmids[i]
            vmax = vmaxs[i]
            cmin = cmins[i]
            cmid = cmids[i]
            cmax = cmaxs[i]
            midnorm = midnorms[i]
            f.write(textwrap.dedent(cptfiletext))
            f.close()
            print('Test loading from CPT file...')
            gmtmap = GMTColorMap.loadFromCPT(tfile)
            print('Passed loading from CPT file...')
            print('Testing vmin/vmax...')
            assert gmtmap.vmin == vmin
            assert gmtmap.vmax == vmax
            print('Passed vmin/vmax...')
            print('Test getting normalized values...')
            norm0 = gmtmap.getNorm(vmin)
            normm = gmtmap.getNorm(vmid)
            norm1 = gmtmap.getNorm(vmax)
            np.testing.assert_almost_equal(norm0,0)
            np.testing.assert_almost_equal(normm,midnorm)
            np.testing.assert_almost_equal(norm1,1)
            print('Passed getting normalized values...')

            print('Testing getting normalized colors at specific values...')
            rgba0 = np.array(gmtmap.getNormColor(vmin))
            rgbam = np.array(gmtmap.getNormColor(vmid))
            rgba1 = np.array(gmtmap.getNormColor(vmax))
            np.testing.assert_almost_equal(rgba0[0,0:3],cmin)
            np.testing.assert_almost_equal(rgbam[0,0:3],cmid)
            np.testing.assert_almost_equal(rgba1[0,0:3],cmax)
            print('Passed getting normalized colors at specific values.')

            print('Testing getting RGB (0-255) colors at specific values...')
            r0,g0,b0,a0 = gmtmap.getRGBColor(vmin)[0]
            rm,gm,bm,am = gmtmap.getRGBColor(vmid)[0]
            r1,g1,b1,a1 = gmtmap.getRGBColor(vmax)[0]
            cmintest = tuple(np.round(np.array(cmin)*255).astype(np.int16))
            cmidtest = tuple(np.round(np.array(cmid)*255).astype(np.int16))
            cmaxtest = tuple(np.round(np.array(cmax)*255).astype(np.int16))
            assert (r0,g0,b0) == cmintest
            assert (rm,gm,bm) == cmidtest
            assert (r1,g1,b1) == cmaxtest
            print('Passed getting RGB (0-255) colors at specific values...')

            # print('Testing save method (and reading back in)...')
            # gmtmap.saveToCPT(tfile)
            # gmtmap2 = GMTColorMap.loadFromCPT(tfile)
            # r0,g0,b0,a0 = gmtmap2.getRGBColor(vmin)[0]
            # rm,gm,bm,am = gmtmap2.getRGBColor(vmid)[0]
            # r1,g1,b1,a1 = gmtmap2.getRGBColor(vmax)[0]
            # cmintest = tuple(np.round(np.array(cmin)*255).astype(np.int16))
            # cmidtest = tuple(np.round(np.array(cmid)*255).astype(np.int16))
            # cmaxtest = tuple(np.round(np.array(cmax)*255).astype(np.int16))
            # assert (r0,g0,b0) == cmintest
            # assert (rm,gm,bm) == cmidtest
            # assert (r1,g1,b1) == cmaxtest
            # print('Passed save method test.')
            if os.path.isfile(tfile):
                os.remove(tfile)
    except Exception as e:
        print('At least one test failed: "%s"' % str(e))
        if os.path.isfile(tfile):
            os.remove(tfile)
    
if __name__ == '__main__':
    test()
