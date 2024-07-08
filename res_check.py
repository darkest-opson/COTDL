# Here the dna seq and list of enzyme to be check for is provided this function will print the map 
import streamlit as st
from Bio.Seq import Seq
from Bio.Restriction import *

def check_res_site_function(ori_seq):
	new_seq = Seq(ori_seq)
	# st.write("Data type of lst:", type(lst))
	rb = RestrictionBatch([], ['C'])
	Ana = Analysis(rb,new_seq, linear=False)
	return Ana
    # res_map = Ana.print_as('map')
    # res_info = Ana.print_that()
    # return res_map, res_info
 