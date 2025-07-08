#!/usr/bin/env python3

# Copyright (C) 2021 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import sys, re, os, json
from argparse import ArgumentParser

mparser = ArgumentParser(usage="Get information on all functionals")
mparser.add_argument('--srcdir', type = str, default = ".",
                     help='Directory where to find the source code')
params = mparser.parse_args()

def read_infos(srcdir, family, all_ids):
  '''Parses a list of files of type family_*.c and hyb_family_*.c,
  reads the list of functionals, ids, and the info structures,
  and returns the information as a dictionary
  '''
  import glob

  # find all files that belong to this family, for example should
  # returns all files of the type gga_*.c for GGAs
  glob_files  = glob.glob(srcdir + "/" + family + "_*.c")
  if family.find('hyb_') != -1:
     glob_files += glob.glob(srcdir + "/" + family.replace('hyb_','') + "_*.c")
  else:
     glob_files += glob.glob(srcdir + "/hyb_" + family + "_*.c")

  # pattern that matches "#define FUNC number /* comment */"
  pattern_def = re.compile(r'#define\s+XC_((' + family.upper() + r")_\S+)\s+(\S+)\s+\/\*\s*(.*?)\s*\*\/")
  # pattern that matches "xc_func_info_type xc_func_info_"
  pattern_inf = re.compile(r'^(const |)xc_func_info_type xc_func_info_(' + family.lower() + r"\S+)")

  numbers  = {}
  infos    = {}

  for mfile in glob_files:
    lines = open(mfile).readlines()
    nline = 0

    while nline < len(lines):
      # match the #define line
      res = pattern_def.match(lines[nline])
      if res is not None:
        name, number = (res.group(1), int(res.group(3)))
        name = name.lower()

        # check if functional number is repeated
        if number in all_ids:
          print("Error: ID '" + number +'" repeated in')
          print("\n" + name + " and in " + all_ids(number))
          sys.exit()

        numbers[name] = number

      # match the xc_func_info_type line
      res = pattern_inf.match(lines[nline])
      if res is not None:
        name = res.group(2)

        struct = ""
        nline += 1
        while nline < len(lines) and lines[nline][0] != "}":
          # remove spaces
          new_line = lines[nline].strip()
          # remove // C comments from line
          new_line = re.sub(r'//.*?', '', new_line)
          # replace commas in braces by |
          new_line = ''.join(m.replace(',', '|') if m.startswith('{') else m
                  for m in re.split(r'(\{[^}]+\})', new_line))
          # replace commas in string by |
          new_line = ''.join(m.replace(',', '|') if m.startswith('"') else m
                  for m in re.split('("[^"]+")', new_line))
          struct += new_line
          nline += 1

        # remove /* ... */ C comments from line
        struct = re.sub(r'/\*.*?\*/', '', struct)

        # split by command not inside braces. This does not work if we have nested braces
        info = [string.replace('|',',') for string in struct.split(",")]

        # build info dictionary
        infos[name] = {}
        infos[name]["number"]   = numbers[name]
        infos[name]["file"]     = os.path.basename(mfile)
        infos[name]["codename"] = name
        infos[name]["kind"]     = info[1]
        infos[name]["description"] = info[2].replace('"', '')
        infos[name]["family"]   = info[3]
        infos[name]["references"] = [i.strip().replace('&xc_ref_', '') for i in re.sub("[{}]", "", info[4]).split('|')]
        infos[name]["dimensionality"] = int(re.sub(r'.*XC_FLAGS_(.)D.*', r'\1', info[5]))
        infos[name]["other_flags"] = [i.strip() for i in info[5].split("|") if not re.search("XC_FLAGS_.D", i)]
        infos[name]["min_dens"] = float(info[6])
        infos[name]["ext_params"] = [i.strip() for i in re.sub("[{}]", "", info[7]).split('|')]
        infos[name]["f_init"]   = info[8]
        infos[name]["f_end"]    = info[9]
        infos[name]["f_lda"]    = info[10]
        infos[name]["f_gga"]    = info[11]
        infos[name]["f_mgga"]   = info[12]
      nline += 1

  return infos

# hybrids are not added automatically
families = ("lda", "hyb_lda", "gga", "hyb_gga", "mgga", "hyb_mgga")
family_infos = {}
all_ids   = {}
all_infos = {}

# Read infos
for family in families:
  infos = read_infos(params.srcdir, family, all_ids)
  family_infos[family] = infos
  all_infos.update(infos)

# check for duplicate IDs
info_list = list(all_infos.values())
for ifun, info in enumerate(info_list):
  for jfun in range(ifun):
    if info_list[ifun]["number"] == info_list[jfun]["number"]:
      number = info_list[ifun]["number"]
      name1 = info_list[ifun]["codename"]
      name2 = info_list[jfun]["codename"]
      print("Error: ID {} repeated in {} and {}".format(number, name1, name2))
      sys.exit()

for family in families:
  full_infos = family_infos[family]
  prefix = ""

  infos = {}
  for func in full_infos:
    infos[func] = full_infos[func]

  # create funcs_family.c file
  fh =  open("funcs_" + prefix + family + ".c", "w")
  fh.write('#include "util.h"\n\n')
  for info in sorted(infos.values(), key=lambda item: item["number"]):
    fh.write('extern xc_func_info_type xc_func_info_' + info["codename"] + ';\n')
  fh.write('\nconst xc_func_info_type *xc_' + prefix + family + '_known_funct[] = {\n')
  for info in sorted(infos.values(), key=lambda item: item["number"]):
    fh.write('  &xc_func_info_' + info["codename"] + ',\n')
  fh.write('  NULL\n};\n')
  fh.close()

# create funcs_key.c file
fh =  open("funcs_key.c", "w")
fh.write('#include "util.h"\n\nxc_functional_key_t xc_functional_keys[] = {\n')
for info in sorted(all_infos.values(), key=lambda item: item["number"]):
  fh.write('  {"' + info["codename"] + '", ' + str(info["number"]) + '},\n')
fh.write('{"", -1}\n};\n')
fh.close()

# create C and F90 files with list of functionals
fh1 =  open("xc_funcs.h", "w")
fh2 =  open("xc_funcs_worker.h", "w")
fh3 =  open("libxc_inc.f90", "w")
for info in sorted(all_infos.values(), key=lambda item: item["number"]):
  line_c = '{:40s} {:5d} {:s}'.format(
    '#define  XC_' + info["codename"].upper(),
    info["number"],
    '/* ' + info["description"] + ' */\n')
  line_f90 = '{:s}\n{:40s} {:s} {:5d}\n\n'.format(
    ' ! ' + re.sub(r'["\']', '', info["description"]),
    ' integer(c_int), parameter, public :: XC_' + info["codename"].upper(),
    " = ",
    info["number"])
  if info["number"] < 100000:
    fh1.write(line_c)
    fh3.write(line_f90)
  else:
    fh2.write(line_c)
fh1.close()
fh2.close()

# dump all information to a json file
with open("libxc_docs.json", "w") as fh:
  json.dump(all_infos, fh, indent=2)
