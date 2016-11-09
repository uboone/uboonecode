<?xml version="1.0" encoding="UTF-8" ?>
<gdml xmlns:gdml="http://cern.ch/2001/Schemas/GDML"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:noNamespaceSchemaLocation="GDMLSchema/gdml.xsd">
  <%include file="define.mako"/>
  <%include file="materials.mako"/>
  <%include file="solids.mako"/>
  <%include file="structure.mako"/>
  <setup name="Default" version="1.0">
    <world ref="volWorld" />
  </setup>
</gdml>
