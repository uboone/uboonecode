#!/usr/bin/python
import copy, os
from mako.template import Template
from mako.lookup import TemplateLookup

from settings import context

if __name__ == "__main__":
    dir_path = os.path.dirname(os.path.realpath(__file__))
    mylookup = TemplateLookup(directories=[os.path.join(dir_path,'templates') ])
    temp = Template('<%include file="microboonevX.mako"/>',
      lookup = mylookup,
      strict_undefined=True)
    with_wires = temp.render(data="world", attributes=context)
    open("microboonevX.gdml",'w').write(with_wires)
