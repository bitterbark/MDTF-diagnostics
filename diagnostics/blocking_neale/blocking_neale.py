# This file is part of the blocking_neale module of the MDTF code package (see LICENSE.txt)
# ============================================================
# Rich Neale's Blocking Code
# Sample code to call NCL from python 
# ============================================================

import os
import sys
import subprocess
import time
import yaml

# ============================================================
# write_cases_page
#
# blocking_neale.html (in POD_HOME) is the static template the framework
# copies to WORK_DIR. The per-case figures depend on how many cases were
# run, so they go on a separate page written straight to WORK_DIR, which
# the framework does not overwrite.
# ============================================================
def write_cases_page(html_page: str, case_dict: dict):
    """write_cases_page: write one row per case, linking its annual cycle figure

    Arguments: html_page (string): file name full path
               case_dict (dict): {case_name: {startdate, enddate, ...}}
    """
    with open(html_page, "w") as f:
        f.write("<html>\n<head>\n<title>Blocking Diagnostic (Rich Neale): MDTF cases</title>\n</head>\n<body>\n")
        f.write('<img src="../mdtf_diag_banner.png">\n')
        f.write("<h3>Blocking Diagnostic (Rich Neale): annual cycle for each MDTF case</h3>\n")
        f.write('<p><A href=blocking_neale.html>Back to blocking_neale results</A></p>\n')
        f.write("<TABLE>\n<TR><TH align=left>Case<TH align=left>Dates<TH align=left>Annual cycle\n")
        for case_name, case_settings in case_dict.items():
            f.write("<TR><TD>{0}<TD>{1} - {2}"
                    "<TD align=center><A href=model/block_freq_anncycle.{0}.png>figure</A>\n".format(
                        case_name, case_settings.get("startdate", ""), case_settings.get("enddate", "")))
        f.write("</TABLE>\n</body>\n</html>\n")


# ============================================================
# generate_ncl_plots - call a nclPlotFile via subprocess call
# ============================================================


def generate_ncl_plots(nclPlotFile):
    """generate_plots_call - call a nclPlotFile via subprocess call
   
    Arguments:
    nclPlotFile (string) - full path to ncl plotting file name
    """
    # check if the nclPlotFile exists - 
    # don't exit if it does not exists just print a warning.
    try:
        pipe = subprocess.Popen(['ncl {0}'.format(nclPlotFile)], shell=True, stdout=subprocess.PIPE)
        output = pipe.communicate()[0].decode()
        print('NCL routine {0} \n {1}'.format(nclPlotFile, output))
        while pipe.poll() is None:
            time.sleep(0.5)
    except OSError as e:
        print('WARNING', e.errno, e.strerror)

    return 0

############################################################
# MAIN
############################################################

# ============================================================
# Translate yaml file variables to environment variables for 
# NCL programs to read
# ============================================================

# Check for $WORK_DIR/case_env.yaml, as sign of multiple cases


print("blocking_neale.py looking for possible multicase case_env_file")
env_var = "case_env_file"
if env_var in os.environ:
    case_env_file = os.environ.get("case_env_file")
    print("blocking_neale.py case_env_file found? ",case_env_file)

    if os.path.isfile(case_env_file):
        with open(case_env_file, 'r') as stream:
            try:
                case_info = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        # case_info.yml holds CASE_LIST and CATALOG_FILE; only CASE_LIST has the cases
        case_info = case_info["CASE_LIST"]

        print("blocking_neale.py read multicase yaml file, setting env vars for NCL")

        icase = 0  # index for cases, needed to save numbered env vars
        for case_name, case_settings in case_info.items():
            icase = icase + 1
            print("case ", icase, ": ", case_name)
            for k, v in case_settings.items():
                casei_env_var_name = "CASE" + str(icase) + "_" + str(k)
                os.environ[casei_env_var_name] = str(v)
                print("setenv ", casei_env_var_name, "\t ", v)

            # NCL reads FIRSTYR/LASTYR; the framework now only provides startdate/enddate
            for yr_var, date_key in (("FIRSTYR", "startdate"), ("LASTYR", "enddate")):
                yr = str(case_settings[date_key])[:4]
                os.environ["CASE" + str(icase) + "_" + yr_var] = yr
                if icase == 1:  # single-case NCL path reads the un-numbered names
                    os.environ[yr_var] = yr
                print("setenv ", "CASE" + str(icase) + "_" + yr_var, "\t ", yr)

        os.environ["CASE_N"] = str(icase)
        print("setenv ", "CASE_N", "\t ", icase)
else:
    print("No multicase case_env_file found so proceeding as single case")
    os.environ["CASE_N"] = "1"

# ============================================================
# Call NCL code here
# ============================================================
if not os.path.exists(os.path.join(os.environ['DATADIR'], 'day')):
    os.makedirs(os.path.join(os.environ['DATADIR'], 'day'))

print("blocking_neale.py calling blocking.ncl")
generate_ncl_plots(os.environ["POD_HOME"]+"/blocking.ncl")

# ============================================================
# Write the per-case page (see write_cases_page)
# ============================================================
if "case_info" not in globals():  # no case_env_file: single case from env vars
    case_info = {os.environ["CASENAME"]: {"startdate": os.environ.get("startdate", ""),
                                          "enddate": os.environ.get("enddate", "")}}
cases_page = os.path.join(os.environ["WORK_DIR"], "blocking_neale_cases.html")
print("blocking_neale.py writing per-case webpage ", cases_page)
write_cases_page(cases_page, case_info)


# ============================================================

print("blocking_neale.py finished.")
sys.exit(0)
