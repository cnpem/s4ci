from PyQt5 import QtWidgets
import create_ui as ui
if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    widget = ui.MyWidget()

    widget.addFormFactors()
    widget.add_dist()
    widget.add_sf()
    widget.addRg()    

    widget.run()
    #ui.sys.exit(app.exec())
    ui.sys.exit(app.exec_())