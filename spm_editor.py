#!/usr/bin/python
"""
Plotting and analysis of CORS station time series

Code based on exmaple
Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 19.01.2009
"""
import sys
import os
import os.path
from datetime import datetime, timedelta
import re
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg
from matplotlib.figure import Figure
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np

import stn_pred_model as spm

help_file='spm_editor_help.html'

model_dir='models'
model_file='{code}_spm.xml'

timeseries_dir='timeseries'
timeseries_file_re=re.compile(r'^(?P<code>\w{4})_igs08_xyz.dat$')
timeseries_file='{code}_igs08_xyz.dat'

class ModelTableView( QTableView ):

    rowSelected = pyqtSignal( int, name="rowSelected" )

    def __init__(self, parent=None):
        QTableView.__init__(self,parent)
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.horizontalHeader().setStretchLastSection(True)
        self.horizontalHeader().setHighlightSections(False)
        self.verticalHeader().setVisible(False)
        self.verticalHeader().setDefaultSectionSize(17)
        # self.setSortingEnabled(True)
        # self.setEditTriggers(QAbstractItemView.AllEditTriggers)
        self.setStyleSheet("* { gridline-color: gray }")

    def selectionChanged( self, selected, deselected ):
        QTableView.selectionChanged( self, selected, deselected )
        self.rowSelected.emit(self.selectedRow())

    def selectedRow( self ):
        rows = self.selectionModel().selectedRows()
        row=-1
        if len(rows) == 1:
            row=rows[0].row()
        return row
    
class ComponentTableModel( QAbstractTableModel ):

    def __init__( self, parent=None ):
        QAbstractTableModel.__init__( self, parent )
        self.model=None
        self._headers=['Use?','Fit?','Name']

    def setPredModel(self, model ):
        self.model=model
        self.modelReset.emit()

    def component(self,row):
        if self.model and row >= 0 and row < len(self.model.components):
             return self.model.components[row]
        return None

    def row( self, component ):
        try:
            row=self.model.components.index(component)
        except:
            row=-1
        return row

    def addComponent(self,component):
        self.beginInsertRows(QModelIndex(), self.rowCount(),self.rowCount())
        self.model.addComponent(component)
        self.endInsertRows()
        self.refresh()

    def removeComponent(self,component):
        try:
            index = self.model.components.index(component)
            self.beginRemoveRows(QModelIndex(),index,index)
            self.model.removeComponent(component)
            self.endRemoveRows()
        except:
            pass

    def refresh(self):
        self.model.sortComponents()
        self.dataChanged.emit(self.index(0,0),self.index(self.rowCount()-1,self.columnCount()-1))

    def count(self):
        return len(self.model.components) if self.model else 0

    def rowCount( self, parent=QModelIndex() ):
        return self.count() if not parent.isValid() else 0

    def columnCount( self, parent=QModelIndex() ):
        return 3 if not parent.isValid() else 0

    def flags( self, index ):
        flag = Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if index.column() < 2:
            flag |= Qt.ItemIsUserCheckable
        return flag

    def data( self, index, role ):
        row = index.row()
        col = index.column()
        component = self.model.components[row]
        if role == Qt.CheckStateRole:
            if col==0:
                return Qt.Checked if component.enabled() else Qt.Unchecked
            if col==1:
                return Qt.Unchecked if component.fixed() else Qt.Checked
        elif role == Qt.DisplayRole or role == Qt.EditRole:
            if col == 2:
                return str(component)
        return QVariant()

    def setData( self, index, value, role ):
        if not index.isValid():
            return False
        row = index.row()
        col = index.column()
        component = self.model.components[row]
        try:
            if role == Qt.CheckStateRole and col==0:
                component.setEnabled(value.toBool())
            elif role == Qt.CheckStateRole and col==1:
                component.setFixed(not value.toBool())
            else:
                return False
            self.dataChanged.emit(index,index)
        except:
            return False
        return True

    def headerData( self, section, orientation, role ):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                if self._headers and section < len(self._headers):
                    return self._headers[section]
        return QVariant()

class ParamTableModel( QAbstractTableModel ):

    def __init__( self, parent=None ):
        QAbstractTableModel.__init__( self, parent )
        self.component=None
        self._headers=['Fit?','Name','Value']

    def setComponent(self, component ):
        self.component=component
        self.modelReset.emit()

    def refresh(self):
        self.dataChanged.emit(self.index(0,0),self.index(self.rowCount()-1,self.columnCount()-1))

    def count(self):
        return len(self.component.parameters) if self.component else 0

    def rowCount( self, parent=QModelIndex() ):
        return self.count() if not parent.isValid() else 0

    def columnCount( self, parent=QModelIndex() ):
        return 3 if not parent.isValid() else 0

    def flags( self, index ):
        flag = Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if index.column() == 2:
            flag |= Qt.ItemIsEditable
        elif index.column() == 0:
            flag |= Qt.ItemIsUserCheckable
        return flag

    def data( self, index, role ):
        row = index.row()
        col = index.column()
        param = self.component.parameters[row]
        if role == Qt.CheckStateRole:
            if col==0:
                return Qt.Unchecked if param.fixed() else Qt.Checked
        elif role == Qt.TextAlignmentRole:
            if col == 2: return Qt.AlignLeft
            if col == 2: return Qt.AlignRight
            if col == 0: return Qt.AlignHCenter
        elif role == Qt.DisplayRole or role == Qt.EditRole:
            if col == 1:
                return param.name()
            elif col == 2:
                return param.getValue()
        return QVariant()

    def setData( self, index, value, role ):
        if not index.isValid():
            return False
        row = index.row()
        col = index.column()
        param = self.component.parameters[row]
        try:
            if role == Qt.EditRole and col==2:
                param.setValue(str(value.toString()))
            elif role == Qt.CheckStateRole and col==0:
                param.setFixed(not value.toBool())
            else:
                return False
            self.dataChanged.emit(index,index)
        except:
            return False
        return True

    def headerData( self, section, orientation, role ):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                if self._headers and section < len(self._headers):
                    return self._headers[section]
        return QVariant()


class NavigationToolbar( NavigationToolbar2QTAgg ):

    def __init__(self, canvas, parent ):
        NavigationToolbar2QTAgg.__init__(self,canvas,parent)
        self.canvas=canvas
        self.clearButtons=[]
        next=None
        for c in self.findChildren(QToolButton):
            if next is None:
                next=c
            if str(c.text()) in ('Subplots','Customize'):
                c.defaultAction().setVisible(False)
                continue
            if str(c.text()) in ('Pan','Zoom'):
                c.toggled.connect(self.clearPicker)
                self.clearButtons.append(c)
                next=None

        pm=QPixmap(32,32)
        pm.fill(QApplication.palette().color(QPalette.Normal,QPalette.Button))
        painter=QPainter(pm)
        painter.fillRect(6,6,20,20,Qt.red)
        painter.fillRect(15,3,3,26,Qt.blue)
        painter.fillRect(3,15,26,3,Qt.blue)
        painter.end()
        icon=QIcon(pm)
        picker=QAction("UseObs",self)
        picker.setIcon(icon)
        picker.setCheckable(True)
        picker.toggled.connect(self.pickerToggled)
        picker.setToolTip("Reject/use observation")
        self.picker = picker
        button=QToolButton(self)
        button.setDefaultAction(self.picker)
        self.insertWidget(next.defaultAction(),button)
    
    def clearPicker( self, checked ):
        if checked:
            self.picker.setChecked(False)

    def pickerToggled( self, checked ):
        if checked:
            for c in self.clearButtons:
                c.defaultAction().setChecked(False)
            self.set_message('Reject/use observation')
            # self.canvas.setCursor(Qt.ArrowCursor)

    def rejectObsMode( self ):
        return self.picker.isChecked()

class AppForm(QMainWindow):

    editableTypes=(spm.equipment_offset,spm.tectonic_offset,spm.slow_slip,spm.slow_slip_ramp,spm.exponential_decay,spm.velocity_change)

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.model=None
        self.undostack=[]
        self.undoptr=-1
        self.autoscale=False
        self.plots=None
        self.obsplots=None
        self.setWindowTitle('PositioNZ station time series analysis')
        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

    def save_plot(self):
        file_choices = "PNG file (*.png)"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def on_about(self):
        msg = """ 
        Station prediction model editor

        View CORS time series and build/edit station prediction models
        """
        QMessageBox.about(self, "About", msg.strip())
    
    def gnsfile(self,code):
        return modelbasefile.replace('{code}',code)

    def checkSaveModel(self):
        model=self.model
        if not model or not model.changed():
            # if model:
            #     print "Model not changed?"
            return
        code=model.station
        result=QMessageBox.question(self,"Save changes to "+code,"Do you want to save changes to "+code+"?",
                                    QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes )
        if result == QMessageBox.Yes:
            self.saveModel()

    def saveModel( self ):
        if self.model:
            self.model.save()

    def savestate( self, clearstack=False ):
        if clearstack:
            self.undostack=[]
        if not self.model:
            return
        xmlstr = self.model.toXmlString()
        if self.undoptr >= 0 and xmlstr == self.undostack[-1]:
            return
        self.undostack.append(xmlstr)
        self.undoptr=len(self.undostack)-1

    def undo( self ):
        if self.model and self.undoptr > 0:
            self.undoptr -= 1
            self.model.loadFromXml( self.undostack[self.undoptr] )
            self.components.refresh()
            self.recalculate(True)

    def redo( self ):
        if self.model and self.undoptr < len(self.undostack)-1:
            self.undoptr += 1
            self.model.loadFromXml( self.undostack[self.undoptr] )
            self.components.refresh()
            self.recalculate(True)

    def reload(self):
        """ Redraws the figure
        """

        code=str(self.codelist.currentItem().text())
        model=self.model
        if model and model.station==code:
            return
        self.checkSaveModel()

        modelFile = self.modelFile(code)
        loadfile=os.path.exists(modelFile)
        self.model=spm.model(station=code,filename=modelFile,loadfile=loadfile)
        self.savestate(True)
        self.model.loadTimeSeries(os.path.join(timeseries_dir,timeseries_file))
        dates,obs,useobs = self.model.getObs()
        message="{0}: {1} observations".format(code,len(dates))
        self.status_text.setText(message)

        self.params.setComponent(None)
        self.eventName.setText('')
        self.components.setPredModel(self.model)
        self.componentTable.resizeColumnsToContents()

        self.calculate()
        self.rescalePlot()
        self.writeStats()

    def loadComponent( self, component ):
        self.removeButton.setEnabled( type(component) in self.editableTypes )
        self.eventName.setText('')
        self.eventName.setEnabled(False)
        self.params.setComponent(component)
        self.paramTable.resizeColumnsToContents()
        if component:
            self.eventName.setText(component.eventName())
            self.eventName.setEnabled(True)
            self.status_text.setText(str(component))

    def writeStats( self ):
        message="Summary stats: {0} days".format(self.obs_count)
        for i,axis in enumerate(['e','n','u']):
            message=message+"\n{0} obs {1:.2f} mm residual {2:.2f} mm".format(
                axis,self.obs_rse[i]*1000,self.residual_rse[i]*1000)
        self.infobox.setPlainText(message)

    def calculate( self ):
        self.obs_rse=[0,0,0]
        self.residual_rse=[0,0,0]
        self.obs_count=0
        dates,obs,useobs = self.model.getObs()
        if dates==None or obs==None or len(obs) == 0:
            return

        self.obs_count=len(dates)
        self.obs_rse=spm.model.robustStandardError(obs)
        residuals=[0,0,0]
        calc=None
        if self.model:
            calc=self.model.calc(dates)
            diff=obs-calc
            # This code was to remove mean difference between observed and calced..
            # Should not be doing this!
            # offset=np.mean(diff,axis=0)
            # calc += offset
            # diff -= offset

            # Calculate a "robust residual" measure based on 95%ile
            residuals=np.percentile(np.abs(diff[useobs]),95.0,axis=0)/(1.96*np.sqrt(2.0))

        self.calc = calc
        self.residual_rse=residuals

    def modelFile(self,code):
        filename=model_dir+'/'+model_file
        filename=filename.replace('{code}',code)
        return filename

    def reloadCodeList(self):
        modelsonly = self.modelsonly.isChecked()
        self.codelist.clear()
        for filename in sorted(os.listdir(timeseries_dir)):
            m=timeseries_file_re.match(filename)
            if not m:
                continue
            code = m.group('code').upper()
            if modelsonly and not os.path.exists(self.modelFile(code)):
                continue
            self.codelist.addItem(code)

    def setAutoscale( self, auto=True ):
        for a in self.plots:
            a.autoscale(auto)

    def rescalePlot( self ):
        self.autoscale=True
        self.replot()

    def recalculate(self, oldstate=False):
        if not oldstate:
            self.savestate()
        self.calculate()
        self.writeStats()
        self.replot()

    def replot(self):
        if not self.model:
            # print "No model loaded!"
            return
        dates,obs,useobs = self.model.getObs()
        if dates==None or obs==None or len(obs) == 0:
            # print "Nothing to plot!"
            return
        plotdays = self.timeasdays.isChecked()
        days=spm.days_array(dates)
        times = days if plotdays else dates


        dmin=np.min(dates)
        dmax=np.max(dates)
        cdates=mdates.num2date(mdates.drange(dmin,dmax+timedelta(0.5),timedelta(1)))
        cdays=spm.days_array(cdates)
        ctimes = cdays if plotdays else cdates

        calc=self.model.calc(cdates)

        plots = self.plots
        for axes in plots:
            if self.autoscale:
                axes.clear()
            else:
                while axes.lines:
                    axes.lines[0].remove()
            if plotdays:
                axes.format_xdata=mticker.ScalarFormatter(useOffset=False)
            else:
                axes.format_xdata=mdates.DateFormatter('%d-%m-%Y')

        title=self.model.station+' time series'
        axis_labels=['East mm','North mm','Up mm']

        self.obsplots=[[None,None],[None,None],[None,None]]
        
        detrend = self.detrend.isChecked()
        trend=0
        for i in range(3):
            if detrend:
                trendp=np.poly1d(np.polyfit(days,obs[:,i],1))
                trend=trendp(days)
            obsi=(obs[:,i]-trend)*1000
            self.obsplots[i][0],=plots[i].plot(times,obsi,'b+',label='Time series',picker=5)
            if detrend:
                trend=trendp(cdays)
            plots[i].plot(ctimes,(calc[:,i]-trend)*1000,marker=None,linestyle='-',color='#00FF00',linewidth=2,label='Prediction model')

            self.obsplots[i][1],=plots[i].plot(times[~useobs],obsi[~useobs],'rs')
            plots[i].set_ylabel(axis_labels[i])
            plots[i].tick_params(labelsize=8)
#        if filename:
#            plt.savefig(filename, bbox_inches='tight')
#        else:
        if True:
            self.figtitle.set_text(title)
            # plt.show()

        if self.autoscale:
            self.setAutoscale(True)
            self.autoscale=False
        if not plotdays:
            self.fig.autofmt_xdate()
        self.canvas.draw()
        self.setAutoscale(False)

    def replotUsed(self):
        if self.obsplots:
            dates,obs,useobs = self.model.getObs()
            for p in self.obsplots:
                dplot=p[0]
                xplot=p[1]
                xplot.set_xdata(dplot.get_xdata()[~useobs])
                xplot.set_ydata(dplot.get_ydata()[~useobs])
        self.canvas.draw()
    
    def updateEventName(self,name):
        component=self.params.component
        if component:
            component.setEventName(str(name))

    def fitComponent(self):
        model=self.model
        component=self.params.component
        if component:
            savefixed = component.fixed()
            component.setFixed(False)
        app=QApplication.instance()
        try:
            app.setOverrideCursor(Qt.WaitCursor)
            ok, message = model.fit()
        finally:
            app.restoreOverrideCursor()
        if component:
            component.setFixed(savefixed)
        self.status_text.setText(message)
        if ok:
            self.components.refresh()
            self.selectComponent(component)
            self.params.refresh()
            self.recalculate()

    def calcDate( self, xvalue ):
        if self.timeasdays.isChecked():
            return spm.fromday(xvalue)
        elif type(xvalue) == datetime:
            return xvalue
        else:
            return mdates.num2date(xvalue)

    def axesClicked( self, event ):
        try:
            self.addComponentDate.setDateTime(self.calcDate(event.xdata))
        except:
            pass

    def obsPicked( self, event ):
        xdata,ydata = event.artist.get_data()
        date=xdata[event.ind][0]
        date=self.calcDate(date)
        self.addComponentDate.setDateTime(date)
        if self.model and self.toolbar.rejectObsMode():
            self.model.setUseObs( event.ind[0], toggle=True )
            self.replotUsed()

    def addComponent( self ):
        ctype = self.addComponentType.itemData(self.addComponentType.currentIndex()).toPyObject()
        cdate = self.addComponentDate.dateTime().toPyDateTime()
        component=ctype(self.model,cdate)
        self.components.addComponent(component)
        self.recalculate()
        self.selectComponent(component)

    def removeComponent( self ):
        component=self.selectedComponent()
        self.components.removeComponent(component)
        self.recalculate()

    def selectedComponent( self ):
        return self.components.component(self.componentTable.selectedRow())

    def selectComponent( self, component ):
        row=self.components.row(component)
        if row >= 0:
            self.componentTable.selectRow(row)

    def create_main_frame(self):
        self.main_frame = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        
        self.dpi = 100
        samescale=False
        title='CORS time series analysis'
        fig, plots=plt.subplots(3,1,sharex=True,sharey=samescale,num=title,figsize=(8,6),dpi=self.dpi)
        self.fig=fig
        self.figtitle=fig.suptitle('No station selected')
        self.plots=plots
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # Other GUI controls

        self.codelist = QListWidget(self)
        self.modelsonly = QCheckBox('Modelled stations only',self)
        self.components = QListWidget(self)

        self.detrend = QCheckBox('Detrend',self)
        self.detrend.setChecked(True)
        self.timeasdays = QCheckBox('Time as days',self)
        self.infobox = QTextEdit(self)
        self.infobox.setReadOnly(True)
        self.components=ComponentTableModel(self)
        self.componentTable=ModelTableView(self)
        self.componentTable.setModel(self.components )
        self.params=ParamTableModel(self)
        self.paramTable=ModelTableView(self)
        self.paramTable.setModel(self.params )
        self.fitButton=QPushButton('Fit',self)
        self.undoButton=QPushButton('Undo',self)
        self.redoButton=QPushButton('Redo',self)
        self.removeButton=QPushButton('Remove',self)
        self.removeButton.setEnabled(False)
        self.addComponentType=QComboBox(self)
        for t in self.editableTypes:
            self.addComponentType.addItem(t.__name__,t)
        self.addComponentDate=QDateTimeEdit(self);
        self.addButton=QPushButton('Add',self)
        self.eventName=QLineEdit(self)
        self.eventName.setEnabled(False)

        self.codelist.itemSelectionChanged.connect( self.reload )
        self.modelsonly.toggled.connect(self.reloadCodeList)
        self.detrend.toggled.connect( self.rescalePlot )
        self.timeasdays.toggled.connect( self.rescalePlot )
        self.componentTable.rowSelected.connect( lambda r: self.loadComponent(self.components.component(r)))
        self.components.dataChanged.connect( lambda x,y: self.recalculate() )
        self.params.dataChanged.connect( lambda x,y: self.recalculate() )
        self.fitButton.clicked.connect( self.fitComponent )
        self.undoButton.clicked.connect( self.undo_action.trigger )
        self.redoButton.clicked.connect( self.redo_action.trigger )
        self.canvas.mpl_connect('button_press_event',self.axesClicked)
        self.canvas.mpl_connect('pick_event',self.obsPicked)
        self.removeButton.clicked.connect( self.removeComponent )
        self.addButton.clicked.connect( self.addComponent )
        self.eventName.textEdited.connect(self.updateEventName)

        #
        # Layout with box sizers
        # 

        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.toolbar,1)
        hbox1.addWidget(self.detrend,0)
        hbox1.addWidget(self.timeasdays,0)

        vbox1 = QVBoxLayout()
        vbox1.addWidget(self.codelist)
        vbox1.addWidget(self.modelsonly)
        vbox1.addWidget(self.infobox)
        
        bbox1 = QHBoxLayout()
        bbox1.addWidget(self.removeButton)
        bbox1.addWidget(self.addComponentType)
        bbox1.addWidget(self.addComponentDate)
        bbox1.addWidget(self.addButton)
        vbox2 = QVBoxLayout()
        vbox2.addWidget(self.componentTable)
        vbox2.addLayout(bbox1)

        nmbox = QHBoxLayout()
        nmbox.addWidget(QLabel('Event',self),0)
        nmbox.addWidget(self.eventName,1)

        bbox2 = QHBoxLayout()
        bbox2.addWidget(self.fitButton)
        bbox2.addWidget(self.undoButton)
        bbox2.addWidget(self.redoButton)

        vbox3 = QVBoxLayout()
        vbox3.addLayout(nmbox)
        vbox3.addWidget(self.paramTable)
        vbox3.addLayout(bbox2)

        hbox2 = QHBoxLayout()
        hbox2.addLayout(vbox1,1)
        hbox2.addLayout(vbox2,2)
        hbox2.addLayout(vbox3,2)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas,2)
        vbox.addLayout(hbox1,0)
        vbox.addLayout(hbox2,1)
        
        self.main_frame.setLayout(vbox)
        self.reloadCodeList()
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("")
        self.statusBar().addWidget(self.status_text, 1)
        
    def show_help( self ):
        helpfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),help_file)
        QDesktopServices.openUrl(QUrl.fromLocalFile(helpfile))

    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        save_action = self.create_action("Save &model",
            shortcut="Ctrl+S", slot=self.saveModel, 
            tip="Save the current model")
        save_plot_action = self.create_action("Save &plot",
            shortcut="Ctrl+P", slot=self.save_plot, 
            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (save_action,save_plot_action, None, quit_action))
        
        self.edit_menu = self.menuBar().addMenu("&Edit")
        self.undo_action = self.create_action("&Undo", slot=self.undo, 
            shortcut="Ctrl+U", tip="Undo the last change")
        self.redo_action = self.create_action("&Redo", slot=self.redo, 
            shortcut="Ctrl+R", tip="Redo the last change")

        self.add_actions( self.edit_menu, (self.undo_action, self.redo_action) )

        self.help_menu = self.menuBar().addMenu("&Help")
        help_action = self.create_action("&Help", 
            shortcut='F1', slot=self.show_help, 
            tip='Brief help information')
        about_action = self.create_action("&About", 
            slot=self.on_about, 
            tip='About the tool')
        
        self.add_actions(self.help_menu, (help_action,about_action))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action

    def closeEvent( self, event ):
        self.checkSaveModel()
        event.accept()


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
