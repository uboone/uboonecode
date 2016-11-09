<volume name="volTPCWireVert">
  <materialref ref="Titanium"/>
  <solidref ref="TPCWireVert"/>
</volume>
<volume name="volTPCPlaneVert">
  <materialref ref="LAr"/>
  <solidref ref="TPCPlaneVert"/>
  % if attributes['wires_on']:
    %for i in range(int( float(attributes['TPCWirePlaneLength']) / float(attributes['TPCWirePitch']) )):
      <physvol>
         <volumeref ref="volTPCWireVert"/>
         <position name="posTPCWireVert${loop.index}" unit="cm"
         z="${attributes['TPCWirePitch']*(-((float(attributes['TPCWirePlaneLength'])/float(attributes['TPCWirePitch'])-1.)/2.+1.)+(i+1.))+0.1}"
         x="0" y="0"/>
         <rotationref ref="rPlus90AboutX"/>
      </physvol>
    %endfor
  %endif
</volume>
