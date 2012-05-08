#
# yanny.py
#
# Python library for reading & writing yanny files.
#
# B. A. Weaver, NYU, 2008-10-20
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
"""
yanny is an object-oriented interface to FTCL/yanny data files following the
specifications at http://www.sdss.org/dr6/dm/flatFiles/yanny.html.  The format
of the returned object is similar to that returned by read_yanny() in the
efftickle perl package (in the photoop product).

Currently multidimensional arrays are only supported for type char.  Also, this
package does not support the <> notation for array size, e.g. float<4>.
"""

__author__ = 'Benjamin Weaver <benjamin dot weaver at nyu dot edu>'

__version__ = '$Revision: 100176 $'.split(': ')[1].split()[0]

__all__ = [ 'yanny', 'read_yanny', 'write_yanny', 'write_yanny_append' ]

#
# Modules
#
import re
import os
import os.path
import datetime
import string

#
# Classes
#
class yanny(dict):
    def __init__(self,*filename):
        """
        Create a yanny object using a yanny file, filename.  If the file exists,
        it is read, & the dict structure of the object will be basically the
        same as that returned by read_yanny() in the efftickle package.
        If the file does not exist, or if no filename is given, a blank
        structure is returned.  Other methods allow for subsequent writing
        to the file.
        """
        #
        # The symbol hash is inherited from the old read_yanny
        #
        self['symbols'] = dict()
        #
        # Create special attributes that contain the internal status of the object
        # this should prevent overlap with keywords in the data files
        #
        self._filename = ''
        self._contents = ''
        #
        # Since the re is expensive, cache the structure types keyed by the field.
        # Create a dictionary for each structure found.
        self.struct_type_caches = dict()
        #
        # If the file exists, read it
        #
        if len(filename) > 0:
            if os.access(filename[0],os.R_OK):
                self._filename = filename[0]
                f = open(filename[0],'r')
                self._contents = f.read()
                f.close()
                self._parse()
        return
    def __str__(self):
        """
        Implement the str function for yanny objects.
        """
        return self._contents
    def __eq__(self,other):
        """
        Test two yanny objects for equality.  The objects are assumed
        to be equal if their contents are equal.
        """
        if isinstance(other,yanny):
            if str(other) == str(self):
                return True
        return False
    def type(self,structure,variable):
        """
        Returns the type of a variable defined in a structure. Returns None
        if the structure or the variable is undefined.
        """
        if structure not in self:
            return None
        if variable not in self[structure]:
            return None
        definition = filter(lambda x: x.find(structure) > 0,
            self['symbols']['struct'])
        if len(definition) != 1:
            return None
		
		# Added code to cache values to speed up parsing large files.
		# 2009.05.11 / Demitri Muna, NYU
        # Find (or create) the cache for this structure.            
        try:
            cache = self.struct_type_caches[structure]
        except KeyError:
            self.struct_type_caches[structure] = dict()
            cache = self.struct_type_caches[structure] # cache for one struct type

        # Lookup (or create) the value for this variable
        try:
            var_type = cache[variable]
        except KeyError:
            typere = re.compile('(\\S+)\\s+%s(\\[.*\\]|);' % variable)
            (typ,array) = typere.search(definition[0]).groups()
            var_type = typ + array
            cache[variable] = var_type

        return var_type
    def isarray(self,structure,variable):
        """
        Returns True if the variable is an array type.  For character types,
        this means a two-dimensional array, e.g.: char[5][20].
        """
        typ = self.type(structure,variable)
        character_array = re.compile(r'char\[.*\]\[.*\]')
        if ((character_array.search(typ) is not None) or
            (typ.find('char') < 0 and typ.find('[') >= 0)):
            return True
        return False
    def convert(self,structure,variable,value):
        """
        Converts value into the appropriate (python) type.  short & int are
        converted to python integer, float & double are converted to python
        float.  Other types are not altered.
        """
        typ = self.type(structure,variable)
        if (typ.find('short') >= 0 or typ.find('int') >= 0):
            if self.isarray(structure,variable):
                return map(int, value)
            else:
                return int(value)
        if (typ.find('float') >= 0 or typ.find('double') >= 0):
            if self.isarray(structure,variable):
                return map(float, value)
            else:
                return float(value)
        return value
    def tables(self):
        """
        Returns a list of all the defined structures.
        """
        foo = self['symbols'].keys()
        foo.remove('struct')
        foo.remove('enum')
        return foo
    def columns(self,table):
        """
        Returns an ordered list of column names associated with a particular
        table.
        """
        foo = list()
        if table in self['symbols']:
            return self['symbols'][table]
        return foo
    def size(self,table):
        """
        Returns the number of rows in a table.
        """
        foo = self.columns(table)
        return len(self[table][foo[0]])
    def pairs(self):
        """
        Returns a list of keys to keyword/value pairs.  Equivalent to doing
        par.keys(), but with all the data tables & other control structures
        stripped out.
        """
        p = list()
        foo = self.tables()
        for k in self.keys():
            if k == 'symbols' or k in foo:
                continue
            p.append(k)
        return p
    def row(self,table,index):
        """
        Returns a list containing a single row from a specified table in
        column order.  If index is out of range, it returns an empty list.
        """
        datarow = list()
        if table in self and index >= 0 and index < self.size(table):
            for c in self.columns(table):
                datarow.append(self[table][c][index])
        return datarow
    def set_filename(self,newfile):
        """
        Updates the filename associated with the yanny object.  Use this if the
        object was created with no filename.
        """
        self._filename = newfile
        return
    def protect(self,x):
        #print type(x)
        if type(x) is str:
            if len(x) == 0 or re.search(r'\s+',x) is not None:
                return '"' + x + '"'
            else:
                return x
        else:
            return str(x)
    def list_of_dicts(self, table):
        """
        Takes a table from the yanny object and constructs a list object
        containing one row per entry. Each item in the list is a dictionary
        keyed by the struct value names.
        """
        return_list = list()
        d = dict()
        
        struct_fields = self.columns(table) # I'm assuming these are in order...
        
        for i in range(self.size(table)):
            one_row = self.row(table, i) # one row as a list
            j = 0
            for key in struct_fields:
                d[key] = one_row[j]
                j = j + 1

            return_list.append(dict(d)) # append a new dict (copy of d)
            
        return return_list
    def new_dict_from_pairs(self):
    	"""
    	Returns a new dictionary (i.e. not a yanny object) based on the keys
    	that self.pairs() returns. There are two reasons this is convenient:
		* the key 'symbol' that is part of the yanny object will not be present
    	* a simple yanny file can be read with no further processing
    		new_dict = yanny(file).new_dict_from_pairs
    	-- added: Demitri Muna, NYU 2009.04.28
    	"""
    	new_dictionary = dict()
    	for key in self.pairs():
    		new_dictionary[key] = self[key]
    	return new_dictionary
    def write(self,*args):
        """
        Write a yanny object to a file.  This assumes that the filename used to
        create the object was not that of a pre-existing file.  If a file of
        the same name is detected, this method will NOT attempt to overwrite
        it, but will print a warning.  This also assumes that the special
        'symbols' key has been properly created.  This will not necessarily
        make the file very human-readable, especially if the data lines are
        long.  If the name of a new file is given, it will write to the new
        file (assuming it doesn't exist).  If the writing is
        successful, the data in the object will be updated.
        """
        if len(args) > 0:
            newfile = args[0]
        else:
            if len(self._filename) > 0:
                newfile = self._filename
            else:
                print "ERROR: No filename specified!"
                return
        basefile = os.path.basename(newfile)
        timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
        contents = "#\n# %s\n#\n# Created by yanny.py\n#\n# %s\n#\n" % (basefile,timestamp)
        #
        # Print any key/value pairs
        #
        for key in self.pairs():
            contents += "%s %s\n" % (key,self[key])
        #
        # Print out enum definitions
        #
        if len(self['symbols']['enum']) > 0:
            contents += "\n" + string.join(self['symbols']['enum'], "\n\n") + "\n"
        #
        # Print out structure definitions
        #
        if len(self['symbols']['struct']) > 0:
            contents += "\n" + string.join(self['symbols']['struct'], "\n\n") + "\n"
        contents += "\n"
        #
        # Print out the data tables
        #
        for sym in self.tables():
            columns = self.columns(sym)
            for k in range(self.size(sym)):
                line = list()
                line.append(sym)
                for col in columns:
                    if self.isarray(sym,col):
                        datum = '{' + string.join(map(self.protect,self[sym][col][k]),' ') + '}'
                    else:
                        datum = self.protect(self[sym][col][k])
                    line.append(datum)
                contents += "%s\n" % string.join(line,' ')
        #
        # Actually write the data to file
        #
        if os.access(newfile,os.F_OK):
            print "%s exists, aborting write!" % newfile
            print "For reference, here's what would have been written:"
            print contents
        else:
            f = open(newfile,'w')
            print >> f, contents
            f.close()
            self._contents = contents
            self._filename = newfile
            self._parse()
        return
    def append(self,datatable):
        """
        Appends data to an existing FTCL/yanny file.  Tries as much as
        possible to preserve the ordering & format of the original file.
        The datatable should adhere to the format of the yanny object, but it
        is not necessary to reproduce the 'symbols' dictionary or the
        '_internal' dictionary.  It will not try to append data to a file that
        does not exist.  If the append is successful, the data in the object
        will be updated.
        """
        if len(self._filename) == 0:
            print "No filename is set for this object. Use the set_filename method to set the filename!"
            return
        if type(datatable) != dict:
            print "Data to append is not of the correct type. Use a dict!"
            return
        timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
        contents = ''
        #
        # Print any key/value pairs
        #
        for key in datatable.keys():
            if key.upper() in self.tables() or key == 'symbols':
                continue
            contents += "%s %s\n" % (key, datatable[key])
        #
        # Print out the data tables
        #
        for sym in self.tables():
            if sym.lower() in datatable:
                datasym = sym.lower()
            else:
                datasym = sym
            if datasym in datatable:
                columns = self.columns(sym)
                for k in range(len(datatable[datasym][columns[0]])):
                    line = list()
                    line.append(sym)
                    for col in columns:
                        if self.isarray(sym,col):
                            datum = '{' + string.join(map(self.protect,datatable[datasym][col][k]),' ') + '}'
                        else:
                            datum = self.protect(datatable[datasym][col][k])
                        line.append(datum)
                    contents += "%s\n" % string.join(line,' ')
        #
        # Actually write the data to file
        #
        if len(contents) > 0:
            contents = ("# Appended by yanny.py at %s.\n" % timestamp) + contents
            if os.access(self._filename,os.W_OK):
                f = open(self._filename,'a')
                print >> f, contents
                f.close()
                self._contents += contents
                self._parse()
            else:
                print "%s does not exist, aborting append!" % self._filename
                print "For reference, here's what would have been written:"
                print contents
        else:
            print "Nothing to be appended!"
        return
    def _parse(self):
        """
        Converts text into tables that users can use.
        """
        def get_token(string):
            """
            Removes the first 'word' from string.  If the 'word' is enclosed
            in double quotes, it returns the contents of the double quotes.
            If the 'word' is enclosed in braces, it returns the contents of the
            braces, but does not attempt to split the array.  If the 'word'
            is the last word of the string, remainder is set equal to the empty
            string.  This is basically a wrapper on some convenient
            regular expressions.
            """
            if string[0] == '"':
                (word, remainder) = re.search(r'^"([^"]*)"\s*(.*)',
                    string).groups()
            elif string[0] == '{':
                (word, remainder) = re.search(r'^\{\s*([^\}]*)\s*\}\s*(.*)',
                    string).groups()
            else:
                try:
                    (word, remainder) = re.split(r'\s+',string,1)
                except ValueError:
                    #print "Problem with string: %s" % string
                    (word, remainder) = (string, '')
            if remainder is None:
                remainder = ''
            return (word,remainder)
        #
        # there are five things we might find
        # 1. 'blank' lines including comments
        # 2. keyword/value pairs (which may have trailing comments)
        # 3. enumeration definitions
        # 4. structure definitions
        # 5. data
        #
        lines = self._contents
        #
        # Reattach lines ending with \
        #
        lines = re.sub(r'\\\s*\n',' ',lines)
        #
        # Find structure & enumeration definitions & strip them out
        #
        self['symbols']['struct'] = re.findall(r'typedef\s*struct\s*\{[^\}]+\}\s*[A-Z_0-9]*\s*;',lines)
        self['symbols']['enum'] = re.findall(r'typedef\s*enum\s*\{[^\}]+\}\s*[A-Z_0-9]*\s*;',lines)
        lines = re.sub(r'typedef\s*struct\s*\{[^\}]+\}\s*[A-Z_0-9]*\s*;','',lines)
        lines = re.sub(r'typedef\s*enum\s*\{[^\}]+\}\s*[A-Z_0-9]*\s*;','',lines)
        #
        # Interpret the structure definitions
        #
        typedefre = re.compile(r'typedef\s*struct\s*\{([^\}]+)\}\s*([A-Z_0-9]*)\s*;')
        for typedef in self['symbols']['struct']:
            typedefm = typedefre.search(typedef)
            (definition,name) = typedefm.groups()
            self[name.upper()] = dict()
            self['symbols'][name.upper()] = list()
            definitions = re.findall(r'\S+\s+\S+;',definition)
            for d in definitions:
                d = d.replace(';','')
                (datatype,column) = re.split(r'\s+',d)
                column = re.sub(r'\[.*\]$','',column)
                self['symbols'][name.upper()].append(column)
                self[name.upper()][column] = list()
        comments = re.compile(r'^\s*#') # Remove lines containing only comments
        blanks = re.compile(r'^\s*$') # Remove lines containing only whitespace
        trailing_comments = re.compile(r'\s*\#.*$') # Remove trailing comments
        if len(lines) > 0:
            for line in lines.split('\n'):
                if len(line) == 0:
                    continue
                if comments.search(line) is not None:
                    continue
                if blanks.search(line) is not None:
                    continue
                #
                # Remove leading & trailing blanks & comments
                #
                line = line.strip()
                line = trailing_comments.sub('',line)
                #
                # Now if the first word on the line does not match a
                # structure definition it is a keyword/value pair
                #
                (key, value) = get_token(line)
                uckey = key.upper()
                if uckey in self['symbols'].keys():
                    #
                    # Structure data
                    #
                    for column in self['symbols'][uckey]:
                        if len(value) > 0 and blanks.search(value) is None:
                            (data,value) = get_token(value)
                            if self.isarray(uckey,column):
                                #
                                # An array value
                                # if it's character data, it won't be
                                # delimited by {} unless it is a multidimensional
                                # string array.  It may or may not be delimited
                                # by double quotes
                                #
                                # Note, we're assuming here that the only
                                # multidimensional arrays are string arrays
                                #
                                arraydata = list()
                                while len(data) > 0:
                                    (token, data) = get_token(data)
                                    arraydata.append(token)
                                self[uckey][column].append(
                                    self.convert(uckey,column,arraydata))
                            else:
                                #
                                # A single value
                                #
                                self[uckey][column].append(
                                    self.convert(uckey,column,data))
                        else:
                            break
                else:
                    #
                    # Keyword/value pair
                    #
                    self[key] = value
        return

#
# Functions
#
def read_yanny(filename):
    """
    Reads the contents of an FTCL/yanny file & returns the data in a hash.
    """
    par = yanny(filename)
    return par.copy()

def write_yanny(filename,datatable):
    """
    Writes the contents of a hash to an FTCL/yanny file.  Ideally used in conjunction with
    read_yanny().
    """
    par = yanny(filename)
    for key in datatable:
        par[key] = datatable[key]
    par.write()
    return

def write_yanny_append(filename,datatable):
    """
    Appends the contents of a hash to an existing FTCL/yanny file.  Ideally used in conjunction with
    read_yanny().
    """
    par = yanny(filename)
    par.append(datatable)
    return

#
# Testing purposes
#
if __name__ == '__main__':
    par = yanny('test.par')
    print par.pairs()
    for p in par.pairs():
        print "%s => %s" % (p, par[p])
    print par.keys()
    print par['symbols'].keys()
    print par['symbols']['struct']
    print par['symbols']['enum']
    print par.tables()
    for t in par.tables():
        print "%s: %d entries" % (t,par.size(t))
        print par.columns(t)
        for c in par.columns(t):
            print "%s: type %s" % (c,par.type(t,c))
            print par[t][c]
    par.write()
    datatable = {'status_update': {'state':['SUCCESS', 'SUCCESS'],
    'timestamp':['2008-06-22 01:27:33','2008-06-22 01:27:36']}}
    par.set_filename('foo.par')
    par.append(datatable)
