import commands,os

for x in xrange(10):

    contents = open('runs%02d.txt' % (x+1),'r').read()

    cmd = 'samweb create-definition prod_good_run_extunbiased_swizzle_v5_inc_dl_v00_p%02d "file_type data and data_tier raw and ub_project.name swizzle_trigger_streams and ub_project.stage mergeext_unbiased and ub_project.version prod_v04_26_04_05 and run_number ' % x
    cmd += contents
    cmd += '"'

    os.system(cmd)


