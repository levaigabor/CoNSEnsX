#!/usr/bin/env python3

from gi.repository import Gtk, Gio
import os

class CsxUI:

    def __init__(self):
        self.builder = Gtk.Builder()
        self.uifile = "csx_gui/csx_gtk.glade"
        self.builder.add_from_file(self.uifile)
        self.window = self.builder.get_object("main_window")

        hb = Gtk.HeaderBar()
        hb.set_show_close_button(True)
        hb.props.title = "CoNSEnsX"
        self.window.set_titlebar(hb)

        button = Gtk.Button()
        icon = Gio.ThemedIcon(name="mail-send-receive-symbolic")
        image = Gtk.Image.new_from_gicon(icon, Gtk.IconSize.BUTTON)
        button.add(image)
        hb.pack_end(button)

        box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        Gtk.StyleContext.add_class(box.get_style_context(), "linked")

        button = Gtk.Button()
        button.add(Gtk.Arrow(Gtk.ArrowType.LEFT, Gtk.ShadowType.NONE))
        box.add(button)

        button = Gtk.Button()
        button.add(Gtk.Arrow(Gtk.ArrowType.RIGHT, Gtk.ShadowType.NONE))
        box.add(button)

        hb.pack_start(box)

        self.window.connect("destroy", lambda q: Gtk.main_quit())

        # INPUT FILES

        self.PDB_FileChooser = self.builder.get_object("PDB_FileChooser")
        self.PDB_idEntry     = self.builder.get_object("PDB_idEntry")
        self.STR_FileChooser = self.builder.get_object("STR_FileChooser")
        self.NOE_FileChooser = self.builder.get_object("NOE_FileChooser")

        # CONTROLS

        self.Karplus1Radio = self.builder.get_object("Karplus1Radio")
        self.Karplus2Radio = self.builder.get_object("Karplus2Radio")
        self.Karplus3Radio = self.builder.get_object("Karplus3Radio")
        self.fitSwitch     = self.builder.get_object("fitSwitch")
        self.FitRangeEntry = self.builder.get_object("FitRangeEntry")
        self.SVD_Switch    = self.builder.get_object("SVD_Switch")
        self.NOE_Switch    = self.builder.get_object("NOE_Switch")
        self.r3AVG_Switch  = self.builder.get_object("r3AVG_Switch")
        self.BicellesRadio = self.builder.get_object("BicellesRadio")
        self.PhagesRadio   = self.builder.get_object("PhagesRadio")

        # BOTTOM BUTTONS

        self.ClearFormButton = self.builder.get_object("ClearFormButton")
        self.StartCsxButton  = self.builder.get_object("StartCsxButton")
        self.WievCSV_Button  = self.builder.get_object("WievCSV_Button")
        self.WievHTML_Button = self.builder.get_object("WievHTML_Button")

        # CONNECT SIGNALS

        self.PDB_FileChooser.connect("file-set", self.on_PDB_file_set)
        self.STR_FileChooser.connect("file-set", self.on_STR_file_set)

        self.Karplus1Radio.connect("toggled", self.on_Karplus_toggled, "1")
        self.Karplus2Radio.connect("toggled", self.on_Karplus_toggled, "2")
        self.Karplus3Radio.connect("toggled", self.on_Karplus_toggled, "3")

        #self.PDB_idEntry.connect("changed", self.on_PDB_entry_edit)

        self.fitSwitch.connect("notify::active", self.on_fitSwitch_activated)

        self.StartCsxButton.connect("clicked", self.on_start_clicked)

        self.BicellesRadio.connect("toggled", self.on_LCmodel_toggled, "bic")
        self.PhagesRadio.connect("toggled", self.on_LCmodel_toggled, "pf1")


        self.window.show_all()


    #def on_PDB_entry_edit(self, widget):
    #    global PDB_id
    #    PDB_id = self.PDB_idEntry.get_text()
    #    print("The selected PDB ID: ", PDB_id)
    #    if PDB_id:
    #        self.PDB_FileChooser.set_sensitive(False)
    #        PDB_file = None
    #    else:
    #        self.PDB_FileChooser.set_sensitive(True)

    #def on_fit_entry_edit(self, widget):
    #    global fit_range
    #    fit_range = self.FitRangeEntry.get_text()
    #    print("The selected fit range: ", fit_range)


    def on_PDB_file_set(self, widget):
        global PDB_file
        PDB_file = self.PDB_FileChooser.get_filename()
        print("The selected PDB file: ", PDB_file)
        self.PDB_idEntry.set_sensitive(False)

    def on_STR_file_set(self, widget):
        global STR_file
        STR_file = self.STR_FileChooser.get_filename()
        print("The selected STR file: ", STR_file)

    #def on_NOE_file_set(self, widget):
    #    global NOE_file
    #    NOE_file = self.NOE_FileChooser.get_filename()
    #    print("The selected NOE file: ", NOE_file)

    karplus_param = 1

    def on_Karplus_toggled(self, widget, number):
        if widget.get_active():
            CsxUI.karplus_param = number
            print("Karplus parameter set: ", CsxUI.karplus_param)

    lc_model = "bic"

    def on_LCmodel_toggled(self, widget, my_lcmodel):
        if widget.get_active():
            CsxUIlcmodel = my_lcmodel

    def on_fitSwitch_activated(self, widget, gparam):
        if widget.get_active():
            self.FitRangeEntry.set_sensitive(True)
        else:
            self.FitRangeEntry.set_sensitive(False)

    def on_start_clicked(self, widget):
        csx_command = "./consensx.py "
        print("./consensx.py ")

        PDB_file = self.PDB_FileChooser.get_filename()
        if PDB_file != None:
            print("-f " + PDB_file)
            csx_command += "-f " + PDB_file + " "
        else:
            print("-p " + self.PDB_idEntry.get_text())
            csx_command += "-p " + self.PDB_idEntry.get_text() + " "

        STR_file = self.STR_FileChooser.get_filename()
        print("-b " + STR_file)
        csx_command += "-b " + STR_file + " "

        NOE_file = self.NOE_FileChooser.get_filename()
        if NOE_file != None:
            print("-r " + self.NOE_FileChooser.get_filename())
            csx_command += "-r " + self.NOE_FileChooser.get_filename() + " "

        print("-d ", CsxUI.karplus_param)
        csx_command += "-d " + str(CsxUI.karplus_param) + " "

        if self.fitSwitch.get_active():
            print("-s ")
            csx_command += "-s "

        fit_range = self.FitRangeEntry.get_text()
        if fit_range != "":
            print("--fit_range " + fit_range)
            csx_command += "--fit_range " + fit_range + " "

        if self.SVD_Switch.get_active():
            print("-R ")
            csx_command += "-R "

        if self.r3AVG_Switch.get_active():
            print("-r3 ")
            csx_command += "-r3 "

        print("-l " + CsxUI.lc_model)
        csx_command += "-l " + CsxUI.lc_model + " "

        print(csx_command)
        os.system(csx_command)

if __name__ == "__main__":
    hwg = CsxUI()
    Gtk.main()
