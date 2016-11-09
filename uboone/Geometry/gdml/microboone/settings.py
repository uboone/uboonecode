from math import pi, sqrt, sin,cos, tan, radians, asin

context = {
  'wires_on':True,
  'cryostat_on':True,
  'pmt_on': True,
  'gen_vetowall':False,
  'gen_crt_a':False,
  'gen_crt_b':True,
  'gen_extras': True,
  'gen_granite': False,
  'DetEnclosureWidth':1483.26,
  'DetEnclosureHeight':1060.,
  'DetEnclosureLength':1483.26,
  'DirtThickness':300.,
  'CryostatInnerRadius':190.5,
  'CryostatOuterRadius':191.61,
  'CryostatLength':1086.49,
  'CryostatEndcapThickness':1.11,
  'CryostatEndcapLength':67.63,
  'UllageLevelFromTop':34.29,
  'TPCDepth':260.,
  'TPCWidth':256.,
  'TPCLength':1045.,
  'TPCWirePitch':0.3,
  'TPCWirePlaneSpacing':0.3,
  'TPCWirePlaneThickness':0.15,
  'TPCWallThickness':1.77,
  'TPCWireThickness':0.015,
  'UVAngle':60.0,
  'NumberOfTPCPlanes':3,
  'InsulationThickness':16.0,
  'RackX':90.5,
  'RackY':212.7,
  'RackZ':60.96,
  'RackThickness':5.8,
  'PMTOffset':0.,
  'UPlaneYOffset':0.,
  'UPlaneZOffset':0.,
  'VPlaneYOffset':0.,
  'VPlaneZOffset':0.,
  'YPlaneZOffset':0.
}

context['SinUVAngle'] = sin( radians(context['UVAngle']) );
context['CosUVAngle'] = cos( radians(context['UVAngle']) );
context['TanUVAngle'] = tan( radians(context['UVAngle']) );

#manual context setting
context['TPCWirePlaneWidth']	=	233.
context['TPCWirePlaneLength']	=	1037.
context['PMTOffset']            =   424.75*2.54/224.75*2.54/224.75*2.54/224.75*2.54/224.75*2.54/224.75*2.54/224.75*2.54/2

context['pi ']  = pi
context['inch'] = 2.54

context['WorldWidth']		=	100.0*context['DetEnclosureWidth']
context['WorldHeight']	=	100.0*context['DetEnclosureHeight']
context['WorldLength']	=	100.0*context['DetEnclosureLength']

context['CathodePlateDepth']	=	0.1
context['CathodeWidth']	=	5.1
context['CathodeHeight']	=	240.0
context['CathodeLength']	=	1042.0

context['GroundPlateWidth']	=	224.0
context['GroundPlateHeight']	=	0.1
context['GroundPlateLength']	=	1100.0

context['tpc_neg_length']	=	1000.0
context['tpc_neg_width']	=	200.0
context['tpc_neg_height']	=	200.0
context['wire_frame_width']	=	9.5
context['wire_plane_height']	=	240.0
context['wire_plane_length']	=	1042.0
context['wires_plength']		=	(context['wire_plane_length'] - 2*context['wire_frame_width'])
context['wires_pwidth']		=	(context['wire_plane_height'] - 2*context['wire_frame_width'])

context['field_cage_width']			=	256
context['field_cage_height']			=	207.6  #increased size to avoid overlap of TubeBottom and Top with TPCActive
context['field_cage_cross_length']	=	sqrt((context['field_cage_width'])**2+(context['field_cage_height']-50)**2)
context['field_cage_length']			=	1011.6 #increased size to avoid overlap of TubeFront and Back with TPCActive
context['field_cage_loop_interval']	=	1 	# =1 is normal, =4 skips 3/4

context['FieldCageTPCClearance']	=	5*context['inch']
context['FieldCageTubeRadius']	=	0.5*context['inch']
context['FieldCageTubeThickness']	=	0.25*context['inch']
context['FieldCageBeamDepth']		=	12.5
context['FieldCageBeamWidth']		=	2.5
context['FieldCageCrossDepth']	=	2.5
context['FieldCageCrossWidth']	=	4
context['FieldCageCrossLength']	=	context['field_cage_cross_length']

context['TPCTotalLength']	=	context['field_cage_length']
context['TPCTotalWidth']	=	context['field_cage_width']
context['TPCTotalHeight']	=	context['field_cage_height']

context['FieldCageLoopLength']	=	context['TPCTotalLength']+2*(context['FieldCageTPCClearance']+2*context['FieldCageTubeRadius'])
context['FieldCageLoopWidth']		=	context['field_cage_width']
context['FieldCageLoopHeight']	=	context['TPCTotalHeight']+2*(context['FieldCageTPCClearance']+2*context['FieldCageTubeRadius'])

context['FieldCageCornerRadius']	=	0.5*context['FieldCageTPCClearance'];
context['FieldCageCornerY']		=	(0.5*context['FieldCageLoopHeight'])-context['FieldCageCornerRadius']-context['FieldCageTubeRadius']
context['FieldCageCornerZ']		=	(0.5*context['FieldCageLoopLength'])-context['FieldCageCornerRadius']-context['FieldCageTubeRadius']

context['FieldCageHeight']	=	context['FieldCageLoopHeight']+2*(context['FieldCageBeamDepth']+context['FieldCageCrossDepth'])
context['FieldCageLength']	=	context['FieldCageLoopLength']+2*(context['FieldCageBeamDepth']+context['FieldCageCrossDepth'])
context['FieldCageWidth']	=	context['field_cage_width']

context['UVWireCount'] = 2400
context['TPCYWirePitch']=context['TPCWirePitch'] / context['SinUVAngle']
context['TPCZWirePitch']=context['TPCWirePitch'] / context['CosUVAngle']

context['NumberWiresPerEdge'] = int( context['TPCWirePlaneWidth'] / context['TPCYWirePitch'] )
context['NumberCenterWires'] = int(context['UVWireCount'] - 2.*context['NumberWiresPerEdge'])

if not context['wires_on']:
    context['NumberWiresPerEdge'] = 0
    context['NumberCenterWires'] = 0

context['aTPC_xos_cathode'] = 0.5*context['CathodeWidth'] + 0.5*context['CathodePlateDepth']
context['aTPC_xos_wires'] = 3.*context['TPCWirePlaneSpacing'] + context['TPCWirePlaneThickness']
context['aTPC_xoffset'] = context['aTPC_xos_cathode'] - context['aTPC_xos_wires']
context['TPCActiveDepth']  = context['TPCDepth'] - context['aTPC_xos_cathode'] - context['aTPC_xos_wires']
context['TPCActiveHeight'] = context['TPCWirePlaneWidth']
context['TPCActiveLength'] = context['TPCWirePlaneLength']-0.2


context['EndcapRmax'] = 0.5*(context['CryostatEndcapLength']+context['CryostatOuterRadius']**2/context['CryostatEndcapLength'])
context['EndcapThetaDeg'] = asin(context['CryostatOuterRadius']/context['EndcapRmax'])*180.0/pi
context['EndcapZcenter'] = 0.5*context['CryostatLength']+context['CryostatEndcapLength']-context['EndcapRmax']
context['EndcapRmin'] = context['EndcapRmax'] - context['CryostatEndcapThickness']
context['UllageThetaDeg'] = asin(context['CryostatInnerRadius']/context['EndcapRmin'])*180.0/pi
context['UllageLength'] = 2*(context['EndcapZcenter'] + sqrt(context['EndcapRmin']**2-context['CryostatInnerRadius']**2))


context['pmt_pos'] = [{'x':-(11.1*2.54+0.005)/2 -141.487, 'y':55.249,'z':108.693 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.342, 'y':55.249,'z':149.287-424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.380, 'y':27.431,'z':72.034 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.387, 'y':-0.303,'z':194.676 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.098, 'y':-28.576,'z':71.407 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.183, 'y':-56.615,'z':108.802 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.2239,'y':-56.203,'z':149.112- 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.363, 'y':54.646,'z':308.909 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.132, 'y':54.693,'z':349.145 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.096, 'y':-0.829,'z':262.947 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -141.031, 'y':-0.706,'z':394.771 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.953, 'y':-56.261,'z':308.572 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.594, 'y':-57.022,'z':349.273 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.929, 'y':55.771,'z':521.066 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.819, 'y':55.822,'z':561.862 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.564, 'y':-0.875,'z':474.028 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.597, 'y':-0.549,'z':606.217 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.540, 'y':-56.323,'z':521.153 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.566, 'y':-56.205,'z':561.549 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.607, 'y':55.800,'z':732.006 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.486, 'y':55.625,'z':772.816 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.570, 'y':-0.051,'z':685.136 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.250, 'y':-0.502,'z':817.141 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.558, 'y':-56.408,'z':732.207 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.550, 'y':-56.284,'z':772.838 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -139.780, 'y':55.822,'z':931.998 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -139.587, 'y':55.313,'z':972.794 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -139.363, 'y':27.607,'z':1010.644 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -140.122, 'y':-0.722,'z':886.531 - 424.75*2.54*.5},
               {'x':-(11.1*2.54+0.005)/2 -139.400, 'y':-28.625,'z':1011.288 - 424.75*2.54*.5},
	           {'x':-(11.1*2.54+0.005)/2-140.004,  'y':-56.309,'z':932.872 - 424.75*2.54*.5},
	           {'x':-(11.1*2.54+0.005)/2-139.721,  'y':-56.514,'z':972.797 - 424.75*2.54*.5}]
context['pmt_paddle_pos']=[{'x':-161.341, 'y':-28.201 + 20/2*2.54,'z':287.161 - 424.75*2.54*.5},
	           {'x':-160.858, 'y':-27.994 + 20/2*2.54,'z':498.501 - 424.75*2.54*.5},
	           {'x':-160.882, 'y':-28.100 + 20/2*2.54,'z':583.333 - 424.75*2.54*.5},
	           {'x':-160.654, 'y':-27.755 + 20/2*2.54,'z':794.575 - 424.75*2.54*.5}]
