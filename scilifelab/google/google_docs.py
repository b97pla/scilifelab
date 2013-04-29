
import gdata.spreadsheet.service
import gdata.docs.service
import gdata
from gdata import MediaSource
from gdata import GDataEntry
import gdata.docs
from scilifelab.google import _from_unicode, _to_unicode

import base64

class GoogleConnection(object):
    
    def __init__(self, credentials):
        self.credentials = credentials       
        self.authenticate(self.credentials)


    def authenticate(self, credentials):

        login, pwd = self._decode_credentials(credentials)
        if not login or not pwd:
            return False

        self.client.email = login
        self.client.password = pwd
        self.client.source = 'SciLifeLab Scripting'
        self.client.ProgrammaticLogin()

        return True

    def _decode_credentials(self, credentials):

        if not credentials:
            return None

        # Split the username and password
        return base64.b64decode(credentials).split(':', 1)

    def get_key(self, object):
        """Get the unique gdocs key identifier for the supplied object"""
        return object.id.text.split('/')[-1].replace('spreadsheet%3A','')

class Document(GoogleConnection):
    
    def __init__(self, credentials):
        # Create a client class which will make HTTP requests with Google Docs server. 
        self.client = gdata.docs.service.DocsService()
        super(Document, self).__init__(credentials)


    def add_spreadsheet(self, name):
        """Create a new spreadsheet with the specified title"""
        new_entry = gdata.GDataEntry()
        new_entry.title = gdata.atom.Title(text=name)
        category = self.client._MakeKindCategory(gdata.docs.service.SPREADSHEET_LABEL)
        new_entry.category.append(category)
    
        ssheet = self.client.Post(new_entry, '/feeds/documents/private/full')
        return ssheet

    def move_to_folder(self, doc, target_folder):
        target_folder = self.get_folder(target_folder)
        if target_folder is None:
            return
        self.client.MoveIntoFolder(doc,target_folder)
        
    def get_folder(self, folder_name):
        """Get a folder if it exists"""
        q = gdata.docs.service.DocumentQuery(categories=['folder'], params={'showfolders': 'true'})
        for entry in (self.client.Query(q.ToUri()).entry or []):
            if _to_unicode(entry.title.text) == _to_unicode(folder_name):
                return entry
    
        return None

class SpreadSheet(GoogleConnection):
    
    def __init__(self, credentials, name=None):
        # Create a client class which will make HTTP requests with Google Docs server. 
        self.client = gdata.spreadsheet.service.SpreadsheetsService()
        super(SpreadSheet, self).__init__(credentials)
        
        self.doc = Document(self.credentials)
        self.ssheet = None
        if name is not None:
            self.ssheet = self.get_spreadsheet(title=name)
            # Create the spreadsheet if it doesn't exist
            if self.ssheet is None:
                self.ssheet = self.doc.add_spreadsheet(name)
    
    def move_to_folder(self, target_folder):
        key = self.get_key(self.ssheet)
        self.doc.move_to_folder(self.ssheet,target_folder)
        self.ssheet = self.get_spreadsheet(key=key)
        
    def get_spreadsheet(self, title=None, key=None):
        """Get an exact match for a spreadsheet"""
        ssheet = None
        if title is not None:
            feed = self.get_spreadsheets_feed(title, True)
            if len(feed.entry) > 0:
                ssheet = feed.entry[0]
        elif key is not None:
            ssheet = self.client.GetSpreadsheetsFeed(key=key)
        return ssheet

    def get_spreadsheets_feed(self, title, exact_match=False):
        """Get a feed of all available spreadsheets, optionally restricted by title.
        """

        # Create a query that restricts the spreadsheet feed to documents
        # having the supplied title.
        q = self._get_query(title, exact_match)
        # Query the server for an Atom feed containing a list of your documents.    
        return self.client.GetSpreadsheetsFeed(query=q)

    def add_worksheet(self, name, rows=0, cols=0, append=False):
        """Add a new worksheet with the specified title to the specified spreadsheet.
        Will overwrite an existing worksheet with the same title unless append is True
        """
        # Check if a worksheet with the same title exists
        ws = self.get_worksheet(name)
        if ws:
            # If we're appending, just return the first object in the feed
            if append:
                return ws
    
            # Otherwise, drop the existing worksheet
            self.client.DeleteWorksheet(ws)
    
        # Add the desired worksheet
        return self.client.AddWorksheet(_to_unicode(name), rows, cols, self.get_key(self.ssheet))
    
    def get_worksheet(self, name):
        """Get an exact match for a worksheet within a spreadsheet"""
        feed = self.get_worksheets_feed(name, True)
        if len(feed.entry) == 0:
            return None
        return feed.entry[0]


    def get_worksheets_feed(self, name=None, exact_match=False):
        """Get a feed of all worksheets in the supplied spreadsheet, \
        optionally restricted by title.
        """
    
        # Create a query that restricts the spreadsheet feed to documents
        # having the supplied title.
        q = self._get_query(name, exact_match)
        # Get the key for the spreadsheet
        k = self.get_key(self.ssheet)
        # Query the server for an Atom feed containing a list of your documents.
        return self.client.GetWorksheetsFeed(key=k, query=q)
    
    
    def write_rows(self, wsheet, header, rows):
        """Write the supplied data rows to the worksheet,
        using the supplied column headers.
        """
        # Get the keys
        ss_key = self.get_key(self.ssheet)
        ws_key = self.get_key(wsheet)
    
        try:
            # As a workaround for the InsertRow bugs with column names,
            # just use single lowercase letters as column headers to start with
            for i in range(0, len(header)):
                self.client.UpdateCell(1, i + 1, chr(97 + i), ss_key, ws_key)
    
            # Iterate over the rows and add the data to the worksheet
            for row in rows:
                row_data = {}
    
                for i, value in enumerate(row):
                    row_data[chr(97 + i)] = unicode(value)
                self.client.InsertRow(row_data, ss_key, ws_key)
    
            # Lastly, substitute the one-letter header for the real string
            for i in range(0, len(header)):
                self.client.UpdateCell(1, i + 1, _to_unicode(header[i]), ss_key, ws_key)
        except:
            return False
    
        return True

    def get_cell_content(self, wsheet, row_start=0, col_start=0, row_end=0, col_end=0):
        """Get the text contents of the cells from the supplied spreadsheet and
        worksheet and from the specified cell range as a two-dimensional list.
        """
    
        if str(row_start) == '0':
            row_start = '1'
        if str(col_start) == '0':
            col_start = '1'
        if str(row_end) == '0':
            row_end = wsheet.row_count.text
        if str(col_end) == '0':
            col_end = wsheet.col_count.text
    
        feed = (self.get_cell_feed(wsheet, row_start, col_start, row_end, col_end) or [])
    
        # Get the dimensions of the 2D-list
        cols = int(col_end) - int(col_start) + 1
        content = []
        for i, cell in enumerate(feed.entry):
            r = i // cols
            c = i - r * cols
            if c == 0:
                row = []
                content.append(row)
            row.append(_to_unicode((cell.content.text or "")))
    
        return content

    def get_cell_feed(self, wsheet, row_start=0, col_start=0, row_end=0, col_end=0):
        """Get a cell feed from the supplied spreadsheet and worksheet and
        from the specified cell range.
        """
    
        if str(row_start) == '0':
            row_start = '1'
        if str(col_start) == '0':
            col_start = '1'
        if str(row_end) == '0':
            row_end = wsheet.row_count.text
        if str(col_end) == '0':
            col_end = wsheet.col_count.text
    
        p = {'min-row': str(row_start),
             'min-col': str(col_start),
             'max-row': str(row_end),
             'max-col': str(col_end),
             'return-empty': 'True'
             }
        query = gdata.spreadsheet.service.CellQuery(params=p)
        return self.client.GetCellsFeed(self.get_key(self.ssheet), self.get_key(wsheet), query=query)
    
    def _get_query(self, title, exact_match):
        """Get a query object for the supplied parameters"""
        
        p = {}
        if title:
            # The urllib.quote method does not handle unicode, so encode the title in utf-8
            p['title'] = _from_unicode(title)
            if exact_match:
                p['title-exact'] = 'true'
            else:
                p['title-exact'] = 'false'
        
        return gdata.spreadsheet.service.DocumentQuery(params=p)
