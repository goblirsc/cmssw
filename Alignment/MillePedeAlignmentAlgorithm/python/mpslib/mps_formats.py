#  Format binary file name by adding prefixes and extension 
#  as defined by one of the supported format descriptions. 
#

from enum import Enum 
import os 

class mps_formats: 
    class mpsFormats(Enum):
        c_compressed = 1 
        c_plain = 2 
        root_local = 3 
        root_xrd = 4 
        plain_text = 5
        unsupported = -1    # identify unsupported file formats

    def __init__(self, eosInstance: str ="eoscms.cern.ch"):
        self.eosInstance = eosInstance
        # map the user-settable formats (key) to the corresponding
        # internal enums 
        self.knownFormats = {
            "Cgz": self.mpsFormats.c_compressed, 
            "C": self.mpsFormats.c_plain, 
            "root": self.mpsFormats.root_local, 
            "xrd": self.mpsFormats.root_xrd, 
            "csv": self.mpsFormats.plain_text, 
        }
        # map the internal enums to the extensions
        # expected by Mille 
        self.suffixes = {
            self.mpsFormats.c_compressed : ".dat.gz",
            self.mpsFormats.c_plain      : ".dat",
            self.mpsFormats.root_local   : ".root",
            self.mpsFormats.root_xrd     : ".root",
            self.mpsFormats.plain_text   : ".csv",
            self.mpsFormats.unsupported  : ".dat",
        }

    # format the base file name stub to the final extension 
    # we expect to find at the end of the Mille job 
    def milleJobTarget(self, fnameBase: str, format: mpsFormats):
        return f"{fnameBase}{self.suffixes[format]}"  
    
    ### format the provided file name stub 
    ### to the format needed in the Mille job 
    def milleOutFile(self, fnameBase: str, format : mpsFormats):
        # if we are targeting compressed C, Mille will 
        # initially write plain C - we  compress afterward
        if format == self.mpsFormats.c_compressed: 
            format = self.mpsFormats.c_plain

        return self.milleJobTarget(fnameBase, format)
    
    # format the provided directory and file name stub 
    # to the path that pede is expected to read
    def pathForPede(self, directory: os.path, fileName: str, format:mpsFormats):
        fileWithExt = self.milleJobTarget(fileName, format) 
        # XRD will read the full path to the file 
        if format == self.mpsFormats.root_xrd:
            return f"root://{self.eosInstance}/"+os.path.join(directory, fileWithExt)
        # all other setups will read a local copy 
        return fileWithExt
        
    # generate the copy command needed to stage the file to the batch node. 
    def copyCommand(self, directory: os.path, fileName: str, format:mpsFormats):
        if format == self.mpsFormats.root_xrd:
            return ""  
        fileWithExt = self.milleJobTarget(fileName, format) 
        return f"xrdcp {os.path.join(directory, fileWithExt)} {fileWithExt}"

    # generate the file compression command if needed
    def compressCommand(self, fileName: str, format: mpsFormats):
        if format == self.mpsFormats.c_compressed:
            return f"gzip {self.milleOutFile(fileName, format)}"
        else:
            return ""
    # generate the file cleanup command if needed
    def cleanupCommand(self, fileName: str, format: mpsFormats):
        if format == self.mpsFormats.root_xrd:
            return ""
        else:
            return f"rm {self.milleJobTarget(fileName, format)}"
    ### Decode a string identifying a 
    ### binary format and return corresponding
    ### Enum. 
    def identifyFormat(self, fmt_str : str):
        try:
            return self.knownFormats[fmt_str] 
        except KeyError:
            print (f"CAREFUL: Binary format '{fmt_str}' is not known to mps_formats.py.\n", 
                   "         Will write plain C binary output. ")
            return self.mpsFormats.unsupported
