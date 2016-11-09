% for i in range(attributes['NumberWiresPerEdge']):
  <volume name="volTPCWire${i}">
    <materialref ref="Titanium"/>
    <solidref ref="TPCWire${i}"/>
  </volume>
% endfor
<volume name="volTPCWireCommon">
  <materialref ref="Titanium"/>
  <solidref ref="TPCWireCommon"/>
</volume>
<volume name="volTPCPlane">
  <materialref ref="LAr"/>
  <solidref ref="TPCPlane"/>
  % for i in range(attributes['NumberWiresPerEdge']):
    <physvol>
      <volumeref ref="volTPCWire${i}"/>
      <position name="posTPCWire${i}" unit="cm"
        y="${-1./4.*(attributes['TPCWirePlaneWidth'] - attributes['TPCWirePlaneLength']/attributes['TanUVAngle']) + attributes['TPCYWirePitch']/2.*(i - (attributes['UVWireCount']-1.)/2.)}"
        z="${1./4.*(attributes['TPCWirePlaneWidth']*abs(attributes['TanUVAngle']) - attributes['TPCWirePlaneLength']) + attributes['TPCZWirePitch']/2.*(i - (attributes['UVWireCount']-1.)/2.)}"
        x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
  % endfor
  % for i in range(attributes['NumberCenterWires']):
  <physvol>
     <volumeref ref="volTPCWireCommon"/>
     <position name="posTPCWire${attributes['NumberWiresPerEdge']+i}"
       unit="cm" y="0"
       z="${attributes['TPCZWirePitch'] * (i - (attributes['NumberCenterWires']-1.)/2. )}"
       x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
  % endfor
  % for i in range(attributes['NumberWiresPerEdge']):
    <physvol>
      <volumeref ref="volTPCWire${attributes['NumberWiresPerEdge'] - i - 1}"/>
      <position name="posTPCWireB${attributes['NumberWiresPerEdge']+attributes['NumberCenterWires']+i}"
        unit="cm"
        y="${1./4.*(attributes['TPCWirePlaneWidth'] - attributes['TPCWirePlaneLength']/attributes['TanUVAngle']) + attributes['TPCYWirePitch']/2*((attributes['NumberWiresPerEdge'] + attributes['NumberCenterWires'] + i-(attributes['NumberWiresPerEdge'] - i - 1.))/2.)}"
        z="${1./4.*(attributes['TPCWirePlaneLength'] - attributes['TPCWirePlaneWidth']*abs(attributes['TanUVAngle'])) + attributes['TPCZWirePitch']/2*((attributes['NumberWiresPerEdge'] + attributes['NumberCenterWires'] + i-(attributes['NumberWiresPerEdge'] - i - 1.))/2.)}"
        x="0"/>
      <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
  % endfor
</volume>
