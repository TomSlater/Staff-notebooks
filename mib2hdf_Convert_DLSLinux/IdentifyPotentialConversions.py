import os
import argparse
import sys

def check_differences(beamline, year, visit, folder = None):
    print("Starting main function")
    # fist check to see which visit is active

    mib_files = []
    raw_dirs = []
    if folder:
        # check that the path folder exists
        raw_location = os.path.join('/dls',beamline,'data', year, visit, 'Merlin', os.path.relpath(folder))
        if not os.path.exists(raw_location):
            print('This folder ', raw_location,'does not exist!') 
            print('The expected format for folder is sample1/dataset1/')
            sys.exit()  # options here?
    else:    
        raw_location = os.path.join('/dls', beamline,'data', year, visit, 'Merlin')
    
    proc_location = os.path.join('/dls', beamline,'data', year, visit, 'processing', 'Merlin')
    if not os.path.exists(proc_location):
        os.mkdir(proc_location)

    # look through all the files in that location and find any mib files
    os.chdir(raw_location)
        
    path_walker = os.walk('.')
    for p, d, files in path_walker:
        # look at the files and see if there are any mib files there
        for f in files:
            if f.endswith('mib'):
                if folder:
                    p = './'+ folder + p[1:]
                mib_files.append((p, f))
                raw_dirs.append(p)

    # print('RAW dirs:  ', raw_dirs)
    # print('MIB files: ', mib_files)

    # look in the processing folder and list all the directories
    proc_dirs = next(os.walk(proc_location))[1]
    # print('Processing dirs:  ', proc_dirs)
    
    # only using the time-stamp section of the paths to compare:
    raw_dirs_check = []
    for folder in raw_dirs:
        raw_dirs_check.append(folder.split('/')[-1]) 
    # compare the directory lists, and see which has the
    to_convert = set(raw_dirs_check) - set(proc_dirs)
    
    print('To CONVERT:  ', to_convert)
    
    return to_convert, mib_files


def main(beamline, year, visit, folder = None):
    check_differences(beamline, year, visit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('beamline', help='Beamline name')
    parser.add_argument('year', help='Year')
    parser.add_argument('visit', help='Session visit code')
    parser.add_argument('folder', nargs= '?', help='Option to add folder')
    v_help = "Display all debug log messages"
    parser.add_argument("-v", "--verbose", help=v_help, action="store_true",
                        default=False)

    args = parser.parse_args()
    
    main(args.beamline, args.year, args.visit, args.folder)