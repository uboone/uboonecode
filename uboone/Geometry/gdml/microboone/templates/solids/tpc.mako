% if attributes['cryostat_on']:
  <box name="TPC"
   lunit="cm"
   x="${attributes['TPCDepth']}"
   y="${attributes['TPCWidth']}"
   z="${attributes['TPCLength']}"/>
  <box name="TPCActive"
   lunit="cm"
   x="${attributes['TPCActiveDepth']}"
   y="${attributes['TPCActiveHeight']}"
   z="${attributes['TPCActiveLength']}"/>
%endif
