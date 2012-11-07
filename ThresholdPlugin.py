# Theshold plugin
#
# Threshold the image.
# Input - one image
# Output - one image
#
# Perry Radau

#History
#Jan. 12, 2003

import vtk
import qt
import os, sys, os.path



# add the plugin directory to the python path
dir,file = os.path.split(__file__)
if dir not in sys.path:
    sys.path.append(dir)

# a QWidget dialog box to show the user interface
class ThresholdPluginInterface(qt.QWidget):
    def __init__(self,parent=None,name=None):
        qt.QWidget.__init__(self,parent,name)
        self.lowerThr = 0
        self.upperThr = 10000
        self.minval = 0
        self.maxval = 10000
        #These synch variables are used to synchronize the slider and spin box updates. Otherwise
        #  a loop occurs where one updates the other which updates the first etc.
        self.synchSliders = False
        self.synchSpinBoxes = False
     
        self.setMinimumSize(50,70)
        self.setGeometry( 100,100,175,100 )
        
        margin = spacing = 5
        self.vbox = qt.QVBoxLayout( self,margin,spacing )

        #Title
        font = qt.QFont("Arial",10)
        font.setBold( True )
        self.lbl0 = qt.QLabel( self, "Label", 0)
        self.lbl0.setText('Threshold Image')
        self.lbl0.setFont( font )
        qt.QToolTip.add( self.lbl0, "Replace pixels with 0 intensity if meeting threshold criteria." )
        self.vbox.addWidget( self.lbl0 )
        #label
        self.lbl = qt.QLabel( self, "Label", 0)
        self.lbl.setText('Lower threshold (intensity)')
        self.vbox.addWidget( self.lbl )
        qt.QToolTip.add( self.lbl, "Pixels with lower intensity will be set to 0." )

        # Create a threshold spin box
        self.sp = qt.QSpinBox( self.minval, self.maxval, 1, self,"Spin Box")
        self.vbox.addWidget( self.sp )
        self.connect( self.sp, qt.SIGNAL("valueChanged(int)"), self.SpinBoxChanged )

        #lower threshold slider
        self.layout0 = qt.QGrid( 2, qt.QGrid.Horizontal, self, "layout" )
        self.layout0.setSpacing( spacing )
        self.sb = qt.QSlider( 0, 100, 1, 0, qt.QSlider.Horizontal, self.layout0, "Slider" )
        self.sliderlbl =  qt.QLabel( self.layout0, "Label", 0)
        self.sliderlbl.setText( "0" )
        self.sliderlbl.setMinimumWidth(28)
        self.AddSlider( self.sb, "Lower threshold (%)",  "Lower threshold as percentage of maximum.", self.sliderlbl )
        self.vbox.addWidget( self.layout0 )
         
        #label
        self.lbl2 = qt.QLabel( self, "Label", 0)
        self.lbl2.setText('Upper threshold (intensity)')
        self.vbox.addWidget( self.lbl2 )
        qt.QToolTip.add( self.lbl2, "Pixels with higher intensity will be set to 0." )

        # Create a threshold spin box
        self.sp2 = qt.QSpinBox( self.minval, self.maxval, 1, self,"Spin Box")
        self.vbox.addWidget( self.sp2 )
        self.connect( self.sp2, qt.SIGNAL("valueChanged(int)"), self.SpinBoxChanged )

        #upper threshold slider
        self.layout = qt.QGrid( 2, qt.QGrid.Horizontal, self, "layout" )
        self.layout.setSpacing( spacing )
        self.sb2 = qt.QSlider( 0, 100, 1, 100, qt.QSlider.Horizontal, self.layout, "Slider" )
        self.sliderlbl2 =  qt.QLabel( self.layout, "Label", 0)
        self.sliderlbl2.setText( "100" )
        self.sliderlbl2.setMinimumWidth(28)
        self.AddSlider( self.sb2, "Upper threshold (%)",  "Upper threshold as percentage of maximum.", self.sliderlbl2 )
        self.vbox.addWidget( self.layout )

    def AddSlider( self, slider, name, tip, sliderlbl ):
        #label
        self.lbl = qt.QLabel( self, "Label", 0)
        self.lbl.setText(name )
        self.vbox.addWidget( self.lbl )
        qt.QToolTip.add( self.lbl, tip )

        slider.setFocusPolicy( 1 )                
        slider.setFixedHeight( self.sb.sizeHint().height() )
        slider.setTickmarks( qt.QSlider.Below )
        slider.setTickInterval( 10 )
        self.connect( slider, qt.SIGNAL("valueChanged(int)"), self.SliderChanged )



    def SpinBoxChanged( self ):
        if self.synchSpinBoxes == True:
            return
        
        #print 'spin box changed '
        self.lowerThr = self.sp.value()
        self.upperThr = self.sp2.value()
         # set the sliders to agree with spin boxes.
        self.synchSliders = True
        low = (self.lowerThr - self.minval)*100.0/(self.maxval - self.minval)
        high = (self.upperThr - self.minval)*100.0/(self.maxval - self.minval)
        self.sb.setValue( low )
        self.sb2.setValue( high )
        self.sliderlbl.setText( str(int(low)) )
        self.sliderlbl2.setText( str(int(high)) )                
        self.synchSliders = False
        self.Calc()
        
    def SliderChanged( self ):
        """ SliderChanged() -- slider has been moved so update thresholds in GUI and then calculate image.
        """
        if self.synchSliders == True:
            return
        
        #print 'slider changed'
        # set the spin boxes to agree with the slider.
        self.lowerThr = self.minval + self.sb.value()*0.01*(self.maxval - self.minval)

        self.upperThr = self.minval + self.sb2.value()*0.01*(self.maxval - self.minval)

        self.synchSpinBoxes = True
        self.sp.setValue( self.lowerThr )
        self.sp2.setValue( self.upperThr )

        #The % labels are set to agree with the slider.
        self.sliderlbl.setText( str( self.sb.value() ) )
        self.sliderlbl2.setText( str( self.sb2.value() ) )
        
        self.synchSpinBoxes = False
        self.Calc()

    def SetRange( self, minval, maxval ):
        self.minval = minval
        self.maxval = maxval

        self.lowerThr = minval + self.sb.value()*0.01*(maxval - minval)
        self.upperThr = minval + self.sb2.value()*0.01*(maxval - minval)

        # set the spin boxes to agree with the slider.
        self.sp.setValue( self.lowerThr )
        self.sp2.setValue( self.upperThr )

        #The % labels are set to agree with the slider.
        self.sliderlbl.setText( str( self.sb.value() ) )
        self.sliderlbl2.setText( str( self.sb2.value() ) )
        
    
    def InitializeGUI( self ):
        """ InitializeGUI() -- set the default values in the GUI.
        """
        #print 'initialize GUI'
        #Determine the maximum intensity of the image.
        try:
            img = self.filter.GetInput()
            range = [0,1]
            img.GetScalarRange( range )
            #print "range ", range
            self.max = range[1]
            self.upperThr = self.max
        except AttributeError:
            pass
        #spin boxes
        self.sp.setValue( self.lowerThr )
        self.sp2.setMaxValue( self.max )
        self.sp2.setValue( self.upperThr )
        self.sp2.setLineStep( 1 )
      
 
    def GetClassName(self):
        return 'ThresholdPluginInterface'

    def Calc( self ):
        """ threshold the image. The result of each application is not cumulative.
        """
        self.filter.SetReplaceIn( False )
        self.filter.SetReplaceOut( True )
        self.filter.ThresholdBetween( self.lowerThr, self.upperThr )
        self.filter.Update()
        #print 'thresholds = ', self.lowerThr, self.upperThr
        
        # call the render function
        if self.filter.renderFunc[0] != None:
            self.filter.renderFunc[0]()



class ThresholdPlugin(vtk.vtkImageThreshold):
    def __init__(self):
        self.SetReplaceOut(0)
        self.SetInValue(0)
        self.ThresholdBetween(0,0)
        self.dialog = ThresholdPluginInterface()
        self.dialog.filter = self
        self.renderFunc = [None]

    def GetNumberOfRequiredInputs(self):
        return 1

    def SetInput(self,input):
        """SetInput()  -- set the input of the plugin.
        """
        # call superclass method
        vtk.vtkImageThreshold.SetInput(self,input)
        
        if input:
            # update the input
            input.UpdateInformation()
            input.SetUpdateExtentToWholeExtent()
            input.Update()
            # get the min and max values
            minval,maxval = input.GetScalarRange()
            minval = min(0,minval)
            maxval = max(1,maxval)
            self.dialog.SetRange(minval,maxval)

    def GetClassName(self):
        """GetClassName()  -- return the name of this plugin
        """
        return "ThresholdPlugin"

    def OCCIShowInterface(self):
        """OCCIShowInterface()  -- display the dialog box

        This method should create the dialog box (if necessary),
        show it, and bring it to the top.
        """
        self.dialog.show()
        self.dialog.raiseW()

    def OCCIHideInterface(self):
        """OCCIHideInterface()  -- hide the dialog box
        """
        self.dialog.hide()

    def OCCIDestroyInterface(self):
        """OCCIDestroyInterface()  -- destroy the interface dialog box
        """
        self.dialog.destroy()
    
    def OCCISetRenderCallback(self,func):
        """OCCISetRenderCallback()  -- set callback to render the image
        """
        self.renderFunc[0] = func








