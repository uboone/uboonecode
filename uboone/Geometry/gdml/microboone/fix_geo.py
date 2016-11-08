#!/usr/bin/python
"""
    Usage: fix_geo.py <filename.gdml>
    consumes filename.gdml and replaces all attributes which evalaute
    numerically with the evaluated expression.
"""

import xml.etree.ElementTree as ET
import logging, ast, sys
import operator as op

# supported operators
operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
             ast.USub: op.neg}

def eval_expr(expr):
    return eval_(ast.parse(expr, mode='eval').body)

def eval_(node):
    if isinstance(node, ast.Num):
        return node.n
    elif isinstance(node, ast.BinOp):
        return operators[type(node.op)](eval_(node.left), eval_(node.right))
    elif isinstance(node, ast.UnaryOp):
        return operators[type(node.op)](eval_(node.operand))
    else:
        raise TypeError(node)

def process_node(node):
    for key, value in node.attrib.iteritems():
        try:
            new_val = eval_expr(value)
            node.attrib[key] = str(new_val)
            if str(value) != str(new_val):
                logging.debug("Found Unevaluated String: "+str(value))
        except Exception:
            continue
    for i in node:
        process_node(i)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logging.debug("Checking File: "+sys.argv[1])
    tree = ET.parse(sys.argv[1])
    root = tree.getroot()
    process_node(root)
    tree.write(sys.argv[1])
