import sys
from PyQt5 import QtCore, QtWidgets, QtGui
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backend_tools import ToolBase, ToolToggleBase
from matplotlib.figure import Figure
from matplotlib.backend_managers import ToolManager
from matplotlib.backends.backend_qt5 import ToolbarQt
import numpy as np
import form_factors as f
import peaks as p
from pathlib import Path
import qtawesome as qta
import structure_factor as s
import lmfit as lmf
import matplotlib.pyplot as plt
import colorsys
plt.rcParams['toolbar'] = 'toolmanager'
from matplotlib import backend_tools
from matplotlib import get_data_path
from re import split
import os

class GroupHideTool(ToolToggleBase):
    """Show lines with a given gid."""
    default_keymap = 'F'
    description = 'Congelar curva atual'
    default_toggled = False
    
    def __init__(self, *args, gid, plot,  **kwargs):
        self.gid = gid
        self.plot = plot

        super().__init__(*args, **kwargs)
        a = qta.icon('mdi.chart-bell-curve-cumulative', options=[{'draw': 'image'}])
        pixmap = a.pixmap(QtCore.QSize(22, 22), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        a = QtGui.QPixmap.toImage(pixmap)

        a.save(get_data_path() + r"\images\freeze_icon.png")
        self.image = get_data_path() + r"\images\freeze_icon.png"
        
    def enable(self, *args):
        self.set_lines_visibility(True)

    def disable(self, *args):
        self.set_lines_visibility(False)

    def set_lines_visibility(self, state):
        for ax in self.figure.get_axes():
            for line in ax.get_lines():
                if line.get_gid() == self.gid:
                    line.set_ydata(self.plot.get_ydata())
                    line.set_xdata(self.plot.get_xdata())
                    line.set_visible(state)
        
        self.figure.canvas.draw()

class XLogScale(ToolToggleBase):
    """X log-scale"""
    default_keymap = 'K'
    description = 'Log-scale on X axes'
    default_toggled = True
    
    def __init__(self, *args,  **kwargs):
        super().__init__(*args, **kwargs)
        
    def enable(self, *args):
        for ax in self.figure.get_axes():
            ax.set_xscale('log')
        self.figure.canvas.draw()

    def disable(self, *args):
        for ax in self.figure.get_axes():
            ax.set_xscale('linear')
        self.figure.canvas.draw()

class YLogScale(ToolToggleBase):
    """X log-scale"""
    default_keymap = 'K'
    description = 'Log-scale on X axes'
    default_toggled = True
    
    def __init__(self, *args,  **kwargs):
        super().__init__(*args, **kwargs)
        
    def enable(self, *args):
        for ax in self.figure.get_axes():
            ax.set_yscale('log')
        self.figure.canvas.draw()

    def disable(self, *args):
        for ax in self.figure.get_axes():
            ax.set_yscale('linear')
        self.figure.canvas.draw()


class MyWidget(QtWidgets.QWidget):
    '''Classe da interface gráfica.'''
    ###################    INITIALIZATION   #########################
    def __init__(self):
        super().__init__()
        #gráfico:
        self._line = [] 
        self.prepare_graph()
        self.ymin_textbox = False

        #parametros simulação:
        self.prepare_params() 
        self.method = []
        self.method.append(0)
        
        #parametros de layout da interface:
        self.layout = QtWidgets.QGridLayout()
        self.left_layout = QtWidgets.QVBoxLayout()
        self.right_layout = QtWidgets.QVBoxLayout()
        self.final_layout = QtWidgets.QHBoxLayout(self)
        self.left_groupbox = QtWidgets.QGroupBox()
        self.right_groupbox = QtWidgets.QGroupBox()
        self.groupbox = QtWidgets.QGroupBox()
        self.gb_layout = []
        self.final_gb_layout = QtWidgets.QVBoxLayout()

        self.adj_tabs = QtWidgets.QTabWidget()
        self.adj_tabs_array = []
        tabs = ['Data', 'Model', 'Fit', 'Peaks']
        for i in range(len(tabs)):
            self.adj_tabs_array.append(QtWidgets.QWidget())
            self.adj_tabs_array[i].layout = QtWidgets.QVBoxLayout()
            self.adj_tabs.addTab(self.adj_tabs_array[i], tabs[i])
            self.prepare_params_for_each_tab(tabs[i])
        self.resize(1700, 950)

        #tabs:
        self.ffs = []
        self.tabs = []
        self.ff_tabs = QtWidgets.QTabWidget()
             
        #arquivos: 
        self.addFile() 
        self.any_file_added = False
        self.first_file_execution = True

        #variáveis de controle
        self.dist_added = False 
        self.sf_added = False        
        self.dist_adj_applied = False
        self.adj = None  
        self.ff_widgets_added = []
        self.dist_widgets_added = []
        self.sf_widgets_added = []  

        #fit
        self.all_fit_record = [] 
        
        self.intensidade_plot = []    

        self.graph_q_unit = 10e-10

    @QtCore.Slot()
    def run(self):
        '''Esta função executa as funções que precisão ser feitas antes da exibição da janela e exibe a jenela.'''
        self.addgraph() #adiciona o gráfico
        #self.adj_graph() #atualiza os valores do grágico e a formatação
        self.left_layout.setAlignment(QtCore.Qt.AlignTop)

        #layout0 = QtWidgets.QVBoxLayout()
        #layout0.setlayout
        layout1 = QtWidgets.QVBoxLayout()
        layout1.addWidget(self.ff_tabs)
        self.adj_tabs_array[0].setLayout(self.file_final_layout)
        self.adj_tabs_array[1].setLayout(layout1)
        self.adj_tabs.setCurrentIndex(1)
        
        #self.final_gb_layout.addWidget(self.ff_tabs)

        #self.groupbox.setLayout(self.final_gb_layout)
        self.left_layout.addWidget(self.adj_tabs)

        self.right_layout.setAlignment(QtCore.Qt.AlignTop)
        self.left_groupbox.setLayout(self.left_layout)
        self.right_groupbox.setLayout(self.right_layout)
        self.left_groupbox.setMaximumWidth(950)
        self.final_layout.addWidget(self.left_groupbox)
        self.final_layout.addWidget(self.right_groupbox)
        self.prepare_fit_widgets()
        self.colors = Colors()
        self.peak_prepare_widgets()
        self.peak_add_tab()
        self.create_lmf_params()
        #lmf.Parameters.update_constraints = self.new_update_constraints
        self.show()    
    
    def prepare_params_for_each_tab(self, tab):
        self.position[tab] = []
        
        self.text_edit_min[tab] = []
        self.text_edit_max[tab] = []
        self.text_edit_value[tab] = []
        self.slider[tab] = []
        self.param_text[tab] = []
        self.param_vary_cb[tab] = []
        self.param_expr_textbox[tab] = []
    
    def prepare_params(self): 
        '''Esta função cria os dicionários onde os sliders e suas caixas de texto serão salvos e a variável self.position'''
        
        #A variável 'self.position' é usada para localização dos sliders no layout. 
        #Cada slider adicionado soma 1 ao valor da variável
        self.position = {}
        
        self.text_edit_min = {}
        self.text_edit_max = {}
        self.text_edit_value = {}
        self.slider = {}
        self.param_text = {}
        self.param_vary_cb = {}
        self.param_expr_textbox = {}
    
    ##################       GRAPH        #######################
    def ylim_changed(self): 
        '''ajusta os limites y do gráfico usando os valor da caixa de texto'''
        min = float(self.ymin_textbox.text()) #recebe o valor mínimo de y
        max = float(self.ymax_textbox.text()) #recebe o valor máximo de y
        self.ax.set_ylim([min, max]) #seta os limites no gráfico
        if len(self._line)>0:
            self._line[-1].figure.canvas.draw() #Plota o grafico atualizado

    def prepare_graph(self): 
        '''Esta função cria as variavéis necessárias para adicionar o gráfico ao layout'''
        self.sim_figure, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.sim_figure) #Cria a fiura
        self.ax.set_title('Scattering Intensity') #define o título do gráfico
        self.ax.set_ylabel('Intensity [u.a.]') #define o texto do eixo y do gráfico
        self.ax.set_xlabel('Q '+'$[\u212B^{-1}]$') #u212B == Angstrom
        self.ax.set_ylim(1e-3,1e12) #limite em y
    
    def addgraph(self): 
        '''Esta função adiciona o gráfico ao layout'''
        self.graph_groupbox = QtWidgets.QGroupBox() #groupbox onde o grafico será adicionado
        self.figure_layout = QtWidgets.QGridLayout() #
        self.tb = ToolManager(self.sim_figure) 
        self.toolbar = ToolbarQt(self.tb, self.graph_groupbox) #barra de ferramentas do grafico

        #self.toolbar.manager.toolbar.add_tool('zoom', 'foo')
        self.figure_layout.addWidget(self.toolbar, 0,0) #barra de ferramentas do gráfico
        self.figure_layout.addWidget(self.canvas,1,0) 
        self.graph_groupbox.setLayout(self.figure_layout)
        self.graph_groupbox.setMinimumHeight(650)
        self.graph_groupbox.setMaximumHeight(650)
        self.right_layout.addWidget(self.graph_groupbox,1)
        backend_tools.add_tools_to_manager(self.tb)
        backend_tools.add_tools_to_container(self.toolbar)

        
        #TODO Essa parte não funciona direito mais. Era para fazer o freeze da curva, mas como agora podem haver mais de uma cruva no gráfico, isso esté errado.
        #self.main_plot, = plt.plot([1, 2, 3], gid='mygroup')
        #self.main_plot.set_visible(False)   
        #self.tb.add_tool('Freeze', GroupHideTool, gid='mygroup', plot = self.main_plot)
        #self.toolbar.add_tool('Freeze', 'foo')

        #Adiciona a opção de mudar a escala do eixo X na barra de ferramentas do gráfico
        self.tb.add_tool('X Log Scale', XLogScale)
        self.toolbar.add_tool('X Log Scale', 'foo')

        #Adiciona a opção de mudar a escala do eixo Y na barra de ferramentas do gráfico
        self.tb.add_tool('Y Log Scale', YLogScale)
        self.toolbar.add_tool('Y Log Scale', 'foo')

        self.ylim() #Adiciona as caixas de texto para mudar os limites do eixo y no gráfico

    def ylim(self):
        '''Esta função cria as caixas de texto para mudar os limites de y no gráfico'''
        ymin, ymax = self.ax.get_ylim() #recebe os limites inferiores e superiores o gráfico
        ymin = ymin+10**np.floor(np.log10(ymin)) 
        ymax = ymax+10**np.floor(np.log10(ymax))

        self.ymax_text = QtWidgets.QLabel('Y maximum: ') 
        self.ymax_textbox = QtWidgets.QLineEdit(self) #Caixa de texto com o valor máximo de y
        self.ymax_textbox.setText(self.format_e(round(ymax, 3))) #Atualiza o texto da caixa
        self.ymax_textbox.setMaximumWidth(80) 
        self.ymax_textbox.setStyleSheet(self.SS()) #Atualiza o texto da caixa
        self.ymax_textbox.returnPressed.connect(self.ylim_changed) #Ao pressionar 'enter' dentro da classe chama a função que ajusta os limites usando os valores das cixas de texto
        self.figure_layout.addWidget(self.ymax_text, 2, 0) 
        self.figure_layout.addWidget(self.ymax_textbox, 2,1) 

        self.ymin_text = QtWidgets.QLabel('Y minimum: ') 
        self.ymin_textbox = QtWidgets.QLineEdit(self) #Caixa de texto com o valor máximo de y
        self.ymin_textbox.setText(self.format_e(ymin)) #Atualiza o texto da caixa
        self.ymin_textbox.setMaximumWidth(80)
        self.ymin_textbox.setStyleSheet(self.SS()) #Atualiza o texto da caixa
        self.ymin_textbox.returnPressed.connect(self.ylim_changed) #Ao pressionar 'enter' dentro da classe chama a função que ajusta os limites usando os valores das cixas de texto
        self.figure_layout.addWidget(self.ymin_text, 3,0)
        self.figure_layout.addWidget(self.ymin_textbox, 3, 1)  

    def updYLim(self): 
        '''Esta função atualiza os valores de y no gráfico baseado na curva plotada e atualiza nas caixas de texto'''
        #TODO preciso rever ela. Não está sendo usada atualmente. Tem alguns problemas em alguns caso.
        ymin, ymax = 10**(np.floor(np.log10(min(self.intensidade_plot)))), 10**(np.floor(np.log10(max(self.intensidade_plot)))+1)
        #print([ymin, ymax])
        self.ax.set_ylim([ymin, ymax])
        if self.ymin_textbox:
            self.ymax_textbox.setText(self.format_e(round(ymax, 3)))
            self.ymin_textbox.setText(self.format_e(ymin))

    def adj_graph(self): 
        '''Esta função ajusta os limites x e y do gráfico''' 
        #TODO Essa função tabém precisa ser revista
        index = self.ff_tabs.currentIndex()
        #self.updYLim()
        xf = (10**np.floor(np.log10(max(self.ffs[index].q)))+0.5)
        xi = 10**(np.floor(np.log10(min(self.ffs[index].q)))+0.5)
        self.ax.set_xlim([xi,xf])
    
    def upd_graph(self): 
        '''Esta função atualiza a curva do grafico'''
        #A variável 'self.intensidade_plot' é o valor do fator de forma simulado multiplicado por uma escala e somado à um background
        self.intensidade_plot = self.sum_intensidade()
        for i in range(len(self._line)):
            self._line[i].set_visible(False)
            for j in self.relative_file:
                if j.currentIndex()-1 == i:
                    self._line[i].set_data(self.ff_file[i].q[self.indmin_s[i]:self.indmax_s[i]], self.intensidade_plot[i][self.indmin_s[i]:self.indmax_s[i]]*self.multiplier_s[i])
                    self._line[i].set_visible(True)
        if len(self._line)>0:
            self._line[-1].figure.canvas.draw()
        self.ax.legend(loc='upper left', prop={'size':8})

    ######################  APPLY DISTRIBUTION AND STRUCTURE FACTOR FUNCTIONS   ################
    def apply_dist(self, index = None, x = [], sigma = None, funcpars = None, vary = True):
        '''Esta função é usada para aplicar a distribuição de tamanhos, caso haja apenas essa opção, em algum dos modelos selecionados na interface quando o valor de algum parêmtro é modificado.
        A variável "index" é usada para indentificar a aba em que a distribuição será aplicada
        A variável "x" é um array com os valores do vetor de espalhamento q
        A variável "sigma" é um float com o valor de sigma da distribuição
        A variável "funcpars" é um dicionário com os parâmetros do modelo
        A variável "vary" é usada para indicar se a intensidade calculada vai ser salva ou não'''

        if index == None: 
            index = self.ff_tabs.currentIndex()
        if len(x)==0:
            x = self.ffs[index].q
        if sigma == None:
            sigma = self.lmf_params['dist'+'_comma_'+'Model'+'_comma_'+str(index)].value
        if funcpars == None:
            funcpars = self.ffs[index].func_pars
        if self.dist_cb[index].isChecked() and sigma!=0:
            if self.dist_box[index].currentIndex() == 0: 
                dist = 'n' 
            if self.dist_box[index].currentIndex() == 1: 
                dist = 'ln' 

            if self.dist_adj_applied:
                intensidade = self.ffs[index].ff_dist(sigma, funcpars, self.param_box[index].currentIndex(), dist, self.min_amp, self.max_amp, self.n)
            else:
                intensidade = self.ffs[index].ff_dist(sigma, funcpars, self.param_box[index].currentIndex(), dist) 
        
        else:
            intensidade = self.ffs[index].ff3(x, funcpars)   
        if vary:
            self.intensidade[index] = intensidade
        else:
            return intensidade

    def apply_sf(self, name, index = False, x = [], funcpars_ff = None, funcpars_sf = None, vary = True):
        '''Esta função é usada para aplicar o fator de estrutura sem considerar uma distribuição de tamanhos em algum dos modelos selecionados na interface, caso haja apenas essa opção,
        quando o valor de algum parâmetro é modificado.
        A variável "name" é um string com o parâmetro do modelo do fator de estrutura que foi modificado
        A variável "index" é usada para indentificar a aba em que a distribuição será aplicada
        A variável "x" é um array com os valores do vetor de espalhamento q
        A variável "funcpars_ff" é um dicionário com os parâmetros do modelo do fator de forma
        A variável "funcpars_sf" é um dicionário com os parâmetros do modelo do fator de estrutura
        A variável "vary" é usada para indicar se a intensidade calculada vai ser salva ou não'''
        if index == None:
            index = self.ff_tabs.currentIndex()
        if funcpars_ff == None:
            funcpars_ff = self.ffs[index].func_pars
        if len(x)==0:
            x = self.ffs[index].q
        if self.sf_cb[index].isChecked():
            if vary:
                self.sf[index].update_all_params(x, funcpars_ff)          
                #self.sf[index].set_Ref(self.sf_ref_expr(index))

                if name in self.sf[index].func_pars.keys():
                    self.lmf_params[name+'_comma_'+'Model'+'_comma_'+str(index)].value = self.slider['Model'][index][name].value()
                    #self.sf[index].func_pars[name]['value'] = self.slider[index][name].value()
                    self.upd_func_pars_sf()
                elif name in self.ffs[index].func_pars.keys():
                    self.ffs[index].ff3(x, self.ffs[index].func_pars)        

                self.sf[index].sf3(x, self.sf[index].func_pars)    

                self.intensidade[index] = self.ffs[index].intensidade*self.sf[index].intensidade
            else:
                int_ff = self.ffs[index].ff3(x, funcpars_ff)        

                int_sf = self.sf[index].sf_funcpars(x, funcpars_sf)    

                return int_ff*int_sf

    def apply_sf_dist(self, name, index = False, x=[], sigma = None, funcpars = None, funcpars_sf = None, vary = True):
        '''Esta função é usada para aplicar o fator de estrutura e uma distribuição de tamanhos em algum dos modelos selecionados na interface caso haja as duas opções
        quando o valor de algum parâmetro é modificado.
        A variável "name" é um string com o parâmetro do modelo do fator de estrutura que foi modificado
        A variável "index" é usada para indentificar a aba em que a distribuição será aplicada
        A variável "x" é um array com os valores do vetor de espalhamento q
        A variável "sigma" é um float com o valor de sigma da distribuição
        A variável "funcpars" é um dicionário com os parâmetros do modelo do fator de forma
        A variável "funcpars_sf" é um dicionário com os parâmetros do modelo do fator de estrutura
        A variável "vary" é usada para indicar se a intensidade calculada vai ser salva ou não'''
        if index == None:
            index = self.ff_tabs.currentIndex()
        if funcpars == None:
            funcpars = self.ffs[index].func_pars
        if funcpars_sf == None:
            funcpars_sf = self.sf[index].func_pars
        if sigma == None:
            sigma = self.lmf_params['dist'+'_comma_'+'Model'+'_comma_'+str(index)].value
        if len(x)==0:
            x = self.ffs[index].q

        if not(self.sf_cb[index].isChecked()) and not(self.dist_cb[index].isChecked()): #não aplica distribuição nem fator de estrutura
            _intensidade_ = self.ffs[index].ff3(x, funcpars) #calcula o fator de forma
            if vary: self.intensidade[index] = _intensidade_
            else: return _intensidade_
        elif not(self.sf_cb[index].isChecked()) and self.dist_cb[index].isChecked(): #aplica apenas distribuição de tamanhos
            self.apply_dist(index, x, sigma, funcpars, vary)
        elif self.sf_cb[index].isChecked() and not(self.dist_cb[index].isChecked()): #aplica apenas fator de estrutura
            self.apply_sf(name, index, x, funcpars, funcpars_sf, vary)
        elif self.sf_cb[index].isChecked() and self.dist_cb[index].isChecked(): 
            if vary: #se a variável 'vary' for True, salva os novos valores nas variáveis correspendentes
                if name in funcpars_sf:
                    self.lmf_params[name+'_comma_'+'Model'+'_comma_'+str(index)].value = self.slider['Model'][index][name].value()
                    #self.sf[index].func_pars[name]['value'] = self.slider[index][name].value()
                    self.upd_func_pars_sf()  #atualiza os valores dos parâmetros dentro do discionario com os valores do fator de estrutura
                self.sf[index].update_all_params(x, funcpars_sf)
            if sigma>0: #verifica se sigma é maior que zero para calcular a distribuição

                #verifica qual distribuição aplicar
                if self.dist_box[index].currentIndex() == 0: 
                    dist = 'n' 
                if self.dist_box[index].currentIndex() == 1: 
                    dist = 'ln' 
                
                #a variável 'self.dist_asj_applied' é booleana e salva se as configuração para o cálculo da distribuição foram alteradas na interface
                #caso tenham sido, passa os valores que podem ter sidos alterados na configuração para a função
                if self.dist_adj_applied:
                    dist_x, dist_y = self.ffs[index].distribuicao(sigma, funcpars, self.param_box[index].currentIndex(), dist, self.min_amp, self.max_amp, self.n)
                else:
                    dist_x, dist_y = self.ffs[index].distribuicao(sigma, funcpars, self.param_box[index].currentIndex(), dist) 

                norm = (np.sum(dist_y)) #normalização

                #verifica qual método será utilizado para aplicar a distribuição de tamanhos e o fator de estrutura no cálculo da intensidade de espalhamento
                if self.sf_aprox_box[index].currentIndex() == 0: #nesse caso será utilizada a 'Decouple Approximation'
                    if vary:
                        self.ffs[index].update_all_params(x, self.ffs[index].func_pars) #atualiza os parâmetros do fator de forma

                    if self.dist_adj_applied:
                        i_dist = self.ffs[index].ff_dist(sigma, funcpars, self.param_box[index].currentIndex(), dist, self.min_amp, self.max_amp, self.n)
                    else:
                        i_dist = self.ffs[index].ff_dist(sigma, funcpars, self.param_box[index].currentIndex(), dist) 

                    dist_x, dist_y = self.ffs[index].dist.return_dist_curve() #calcula o 'x' e 'y' da distribuição

                    int_dist = np.sum(dist_y)/norm 

                    int_sqrt = [] 
                    for i in dist_x: #calculo do fator de forma para cada valor da distribuição aplicada
                        int_sqrt.append(np.sqrt(self.ffs[index].ff2(i, funcpars, self.dist_box[index].currentIndex()))) #calculo da intensidade

                    int_sqrt = np.transpose(int_sqrt)
                    a_ = (int_sqrt*dist_y)/norm #integrando
                    i_ = np.sum(a_,axis=1) #calculo da integral
                    
                    #if vary: #self.sf[index].set_Ref(self.sf_ref_expr(index))
                    
                    _intensidade_ = i_dist + (1/int_dist)*i_**2*(self.sf[index].sf_funcpars(x, funcpars_sf))

                elif self.sf_aprox_box[index].currentIndex() == 1 and self.param_box[index].currentIndex()>=0: #nesse caso é usado o método 'Local Monodisperse Approximation'
                    if vary:
                        self.ffs[index].update_all_params(x, self.ffs[index].func_pars) #atualiza os valores dos parâmetros dentro do discionario com os valores do fator de estrutura
                    norm = (np.sum(dist_y)) #normalização
                    intensidade = [] 
                    r_ef = self.slider['Model'][index]['r_ef']  #valor de raio efetivo
                    #print('r_ef: {}, ff.p: {}, expr: {}, param: {}'.format(r_ef, self.ff.p, self.expr, self.param_box.currentIndex()))
                    if self.ffs[index].p in self.expr[index] and self.param_box[index].currentIndex()>0:
                        for i in range(len(dist_x)): #calculo do fator de forma para cada valor da distribuição aplicada
                            _fp_sf_ = funcpars_sf
                            _fp_sf_['r_ef']['value'] = r_ef[i]
                            sf_ = self.sf[index].sf_funcpars(x, _fp_sf_)
                            #print(self.ff.current_method)
                            #print(self.methodbox.currentIndex())
                            i_ = self.ffs[index].func(dist_x[i], self.ffs[index].current_method, False)
                            intensidade.append(self.ffs[index].ff2(i_, funcpars, self.param_box[index].currentIndex())*sf_) #calculo da intensidade
                            print('x: {}, dist_x: {}, r_ef:{}, param:{}'.format(i_,dist_x[i], r_ef[i], self.param_box[index].currentIndex()))
                    elif self.ffs[index].p in self.expr[index] and self.param_box[index].currentIndex()==0:
                        for i in range(len(dist_x)):
                            _fp_sf_ = funcpars_sf
                            _fp_sf_['r_ef']['value'] = r_ef[i]
                            sf_ = self.sf[index].sf_funcpars(x, _fp_sf_)
                            #print(self.ff.current_method)
                            #print(self.methodbox.currentIndex())
                            intensidade.append(self.ffs[index].ff2(dist_x[i], funcpars, self.param_box[index].currentIndex())*sf_) #calculo da intensidade
                            print('x: {}, dist_x: {}, r_ef:{}, param:{}'.format(dist_x[i],dist_x[i], r_ef[i], self.param_box[index].currentIndex()))
                    else:
                        for i in range(len(dist_x)):
                            intensidade.append(self.ffs[index].ff2(dist_x[i], funcpars, self.param_box[index].currentIndex())*self.sf[index].sf_funcpars(self.ffs[index].q, funcpars_sf)) #calculo da intensidade

                    intensidade = np.transpose(intensidade)
                    a_ = (intensidade*dist_y)/norm #integrando
                    i_ = np.sum(a_, axis=1) #calculo da integral
                    #print(i_)
                    _intensidade_ = i_

                elif self.sf_aprox_box[index].currentIndex == 2:
                    pass
                    
                if vary:
                    self.intensidade[index] = _intensidade_
                else: return _intensidade_

            else:
                self.apply_sf(name, index, x, funcpars, funcpars_sf, vary)         
    
    def upd_func_pars_ff(self):
        '''Esta função atualiza os valores de todos dos parâmetros de todos os dicionários de modelos de fatores de forma '''
        for index in range(len(self.ffs)):
            for name in self.ffs[index].func_pars:
                #print(name+str(index), self.lmf_params[name+str(index)].value)
                self.ffs[index].func_pars[name]['value'] = self.lmf_params[name+'_comma_'+'Model'+'_comma_'+str(index)].value
                #print(name+str(index), ' :',self.ffs[index].func_pars[name]['value'])

    def upd_func_pars_sf(self):
        '''Esta função atualiza os valores de todos dos parâmetros de todos os dicionários de modelos de fatores de estrutura'''
        #print('upd_fun_pars_sf')
        for index in range(len(self.ffs)):
            if self.sf_cb[index].isChecked():
                for name in self.sf[index].func_pars:
                    self.sf[index].func_pars[name]['value'] = self.lmf_params[name+'_comma_'+'Model'+'_comma_'+str(index)].value
    
    def funcionalities_validation(self, index):
        var, var1, var2, var3 = False, False, False, False
        #TODO revisar essa parte depois, por enquanto funciona mas não sei se está correto para todos o casos
        if self.ff_widgets_added[index]:
            var1 = True
        if self.dist_added:
            if self.dist_widgets_added[index]:
                var2 = True
        if self.sf_added:
            if self.sf_widgets_added[index]:
                var3 = True
        if var1 and var2 and var3:
            var = True
        return var
    
    def calculate(self, name, index, upd, tab):
        '''Essa função recalcula o valor da intensidade de espalahamento de um dos modelos.
        A variável "name" é uma string informanda qual parâmetro do modelo mudou de valor
        A variável "index" é usada para indicar o modelo referente a qual aba deve ser recalculado
        A variável "upd" é booleana e informa se o gráfico da interface deve ser atualizado ou não.
        A variável "tab" informa à qual aba o parâmetro que foi variado pertence ('Model' ou 'Peaks')'''
        self.upd_func_pars_ff() #atualiza os valores dos parâmetros do fator de forma
        self.upd_func_pars_sf() #atualiza os valores dos parâmetros do fator de estrutura
        if tab == 'Peaks' and self.peaks_activate_cb[index].isChecked(): #nesse caso a aba a qual o parâmetro variádo pertence é a aba 'Peaks'
            self.upd_peaks(index) #atualiza a posição dos picos simulados no gráfico
        elif not(self.dist_added) and not(self.sf_added) and self.ff_widgets_added[index]:
            self.intensidade[index] = self.ffs[index].ff3(self.ffs[index].q, self.ffs[index].func_pars) #se não for aplicada distribuição ou fator de estrutura, apenas calcula o fator de forma normalmente
        elif self.dist_added and not(self.sf_added) and self.ff_widgets_added[index] and self.dist_widgets_added[index]: #aplica distribuição
            self.apply_dist(index)
        elif not(self.dist_added) and self.sf_added and self.ff_widgets_added[index] and self.sf_widgets_added[index]: #aplica fator de estrutura
            self.apply_sf(name, index)
        elif self.dist_added and self.sf_added and self.ff_widgets_added[index] and self.dist_widgets_added[index] and self.sf_widgets_added[index]: #aplica distribuição e fator de estrutra
            self.apply_sf_dist(name, index)                
            
        if upd:
            self.upd_graph() #atualiza a curva do grafico
    
    def valuechange(self, name, index = None, upd = True, tab = 'Model', calc = True): 
        '''Esta  é função chamada quando o valor de qualquer slider muda.
        A variável "name" é uma string informanda qual parâmetro do modelo mudou de valor
        A variável "index" é usada para indicar o modelo referente a qual aba deve ser recalculado
        A variável "upd" é booleana e informa se o gráfico da interface deve ser atualizado ou não.
        A variável "tab" informa à qual aba o parâmetro que foi variado pertence ('Model' ou 'Peaks')
        A variável "calc" é booleana e informa se a intensidade de espalhamente deve ser recalculada ou não'''
        if index == None:
            index = self.ff_tabs.currentIndex()

        var = self.funcionalities_validation(index)

        if var:
            # Em todos os casos aqui, a consição: 'if self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].expr == None:' verifica se o atributo 'expr' do parâmetro e nulo para salvar 
            #o valor do slider no parâmetro do lmfit
            if name == 'scale': #verifica se o slider referente à escala foi alterado
                #self.scale[index] = 10**self.slider[index][name].value() #salva o valor na variável referente à escala do fator de forma
                if self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].expr == None:
                    self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].value = self.slider[tab][index][name].value()
            elif name == 'dist':
                if self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].expr == None:
                    self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].value = self.slider[tab][index][name].value()
            elif name == 'peak_q':
                if self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].expr == None:
                    self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].value = self.slider[tab][index][name].value()
            elif name in self.ffs[index].func_pars.keys():
                #verifica se o slider alterado é referente à algum dos dos parâmetros do fator de forma
                #salva o valor do slider modificado na variável func_pars do objeto fator de forma
                if self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].expr == None:
                    self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].value = self.slider[tab][index][name].func_value(self.slider[tab][index][name].value(), self.method[index])
                #requires_update = {name for name, par in self.lmf_params.items() if par._expr is not None}
                #print(requires_update, self.lmf_params[name+str(index)]._expr_deps)
            self.lmf_params._asteval.symtable[name+str(index)] = self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].value #faz com que seja possível usar a string 'name+str(index)' para relacionar as vairáveis atrávés do lmft expressions
            self.lmf_params.update_constraints() #atualiza todos os valores dos parâmetros do lmfit levando em as dependências entre si
                
            
            '''for i in range(len(self.tabs)):
                print('Tab: {}, ff_w: {}, sf_w: {}, dist_w: {}'.format(i, self.ff_widgets_added[i], self.sf_widgets_added[i], self.dist_widgets_added[i]))'''
            
            if calc:
                self.calculate(name, index, upd, tab)

            #verifica se o valor nas caixas de texto dos sliders tem módulo menor que 1000. Caso maior aplica notação cientifica
            if name in self.slider[tab][index].keys():
                if abs(self.slider[tab][index][name].value()) < 1000: 
                    self.text_edit_value[tab][index][name].setText(str(self.slider[tab][index][name].value()))
                else: 
                    self.text_edit_value[tab][index][name].setText(self.format_e(self.slider[tab][index][name].value()))      

    def update_slider(self, name, _min, _max, _value ): 
        '''Esta função atualiza o valor mínimo máximo e atual dos slider com a chave salva na variável 'name' baseado nos valores da caixa de texto'''
        index = self.ff_tabs.currentIndex()
        tab = self.adj_tabs.tabText(self.adj_tabs.currentIndex())
        self.slider[tab][index][name].setMinimum(_min)
        self.slider[tab][index][name].setMaximum(_max)
        self.slider[tab][index][name].setValue(_value)
        self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].set(min=_min)
        self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].set(max=_max)
        self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].set(value=_value)

    #####################   GENERAL FUNCTIONS     ###########################    
    def addparam(self, name, min, max, init_value, symbol, layout = QtWidgets.QGridLayout(), precision = 3, tab = 'Model'):  
        '''Esta função adicona um slider
        A variável 'name' serve como chave no dicionário de sliders
        A variável 'min' define o valor mínimo do slider
        A variável 'max' define o valor máximo do slider
        A variável 'init_value' define o valor inicial do slider'''

        if tab == 'Model':
            index = self.ff_tabs.currentIndex()
        elif tab == 'Peaks':
            index = self.peaks_tabs.currentIndex()

        #print(len(self.param_text)-1, self.ff_tabs.currentIndex())
        self.param_text[tab][index].update({name:QtWidgets.QLabel(symbol)}) #dicionário de labels 
        self.param_vary_cb[tab][index].update({name:QtWidgets.QCheckBox()}) #dicionário de checkboxs
        if name == 'scale':
            self.param_vary_cb[tab][index][name].setChecked(False)
        else:
            self.param_vary_cb[tab][index][name].setChecked(True)

        self.param_expr_textbox[tab][index].update({name:QtWidgets.QLineEdit()})
        #self.text_edit_min[index][name].setMaximumWidth(80)

        #cria as caixas de texto para os valores máximos, mínimos e atuais
        self.text_edit_min[tab][index].update({name:QtWidgets.QLineEdit(str(min))}) #dicionário das caixas de texto dos valores mínimos
        self.text_edit_min[tab][index][name].setMaximumWidth(80)
        
        self.text_edit_max[tab][index].update({name:QtWidgets.QLineEdit(str(max))}) #dicionário das caixas de texto dos valores máximos
        self.text_edit_max[tab][index][name].setMaximumWidth(80)

        self.text_edit_value[tab][index].update({name:QtWidgets.QLineEdit(str(init_value))}) #dicionário das caixas de texto dos valores atuais
        self.text_edit_value[tab][index][name].setMaximumWidth(100)
        
        #adiciona o slider no dicionario de sliders
        self.slider[tab][index].update({name:DoubleSlider(precision, QtCore.Qt.Horizontal)}) #dicionário de sliders
        self.slider[tab][index][name].setMinimum(min)
        self.slider[tab][index][name].setMaximum(max)
        self.slider[tab][index][name].setValue(init_value)
        self.slider[tab][index][name].setMinimumWidth(400)
        self.slider[tab][index][name].setMaximumWidth(600)

        self.toolTip(name, tab) #adiciona legenda ao passar o mouse por cima
        self.alignment(name,'center', tab) #alinhamento do texto

        #adiciona os elementos ao layout
        layout.addWidget(self.param_vary_cb[tab][index][name], self.position[tab][index], 5)
        layout.addWidget(self.param_text[tab][index][name], self.position[tab][index], 6)
        layout.addWidget(self.slider[tab][index][name], self.position[tab][index], 7)
        layout.addWidget(self.text_edit_min[tab][index][name], self.position[tab][index], 8)
        layout.addWidget(self.text_edit_max[tab][index][name], self.position[tab][index], 9)
        layout.addWidget(self.text_edit_value[tab][index][name], self.position[tab][index], 10)
        layout.addWidget(self.param_expr_textbox[tab][index][name], self.position[tab][index], 11)   

        #funções a serem chamadas quando os valores das ciaxas de texto são alterados  
        self.text_edit_min[tab][index][name].returnPressed.connect(lambda: self.update_slider(name, float(self.text_edit_min[tab][index][name].text()),
                                                                                    float(self.text_edit_max[tab][index][name].text()),
                                                                                    float(self.text_edit_value[tab][index][name].text())))
        self.text_edit_max[tab][index][name].returnPressed.connect(lambda: self.update_slider(name, float(self.text_edit_min[tab][index][name].text()),
                                                                                    float(self.text_edit_max[tab][index][name].text()),
                                                                                    float(self.text_edit_value[tab][index][name].text()) ))
        self.text_edit_value[tab][index][name].returnPressed.connect(lambda: self.update_slider(name, float(self.text_edit_min[tab][index][name].text()),
                                                                                    float(self.text_edit_max[tab][index][name].text()),
                                                                                    float(self.text_edit_value[tab][index][name].text()) ))
        self.param_expr_textbox[tab][index][name].returnPressed.connect(lambda: self.param_expr_changed(name, self.ff_widgets_added))
        self.slider[tab][index][name].valueChanged.connect(lambda: self.valuechange(name, None, True, tab)) #função chamada quando o valor do slider muda

        self.position[tab][index]+=1 #varivael para formatação visual

        #formatação visual
        self.text_edit_value[tab][index][name].setStyleSheet(self.SS())
        self.text_edit_min[tab][index][name].setStyleSheet(self.SS())
        self.text_edit_max[tab][index][name].setStyleSheet(self.SS())
        self.param_expr_textbox[tab][index][name].setStyleSheet(self.SS())
        
        return layout

    def param_expr_changed(self, name, finished):
        '''Essa função e usada para trocar a expressão de um determidado parâmetro do modelo.
        A variável "name" é uma string informanda qual parâmetro do modelo mudou de valor
        A variável "finished" é booleana e informa se as alterações na interface gráfica enquanto o programa está sendo executado ja foram finalizadas.
        Ela é necessária pois por conta dessas alterações essa função pode ser chamada em momentos que ela não deveria.'''
        index = self.ff_tabs.currentIndex()
        tab = self.adj_tabs.tabText(self.adj_tabs.currentIndex())
        #print(self.lmf_params[name+str(index)].vary)
        #print(self.lmf_params.valuesdict())
        if self.param_expr_textbox[tab][index][name].text() == '':
            self.slider[tab][index][name].setEnabled(True)
            if finished:
                self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].set(expr='')
        else:
            self.slider[tab][index][name].setEnabled(False)
            if finished:
                self.lmf_params[name+'_comma_'+tab+'_comma_'+str(index)].set(expr=self.param_expr_textbox[tab][index][name].text())
        self.valuechange(name)
        #print(self.lmf_params[name+str(index)].vary)        
    
    def addNeededWidgets(self, tab = 'Model'):
        '''Essa função cria os dicinários dos widgets para adicionar uma nova aba na interface.
        A variável "tab" informa à qual aba o parâmetro que foi variado pertence ('Model' ou 'Peaks')'''
        self.position[tab].append(0)

        self.text_edit_min[tab].append({})
        self.text_edit_max[tab].append({})
        self.text_edit_value[tab].append({})
        self.slider[tab].append({})
        self.param_text[tab].append({})
        self.param_vary_cb[tab].append({})
        self.param_expr_textbox[tab].append({})
    
    def sum_intensidade(self):
        '''Essa função soma a intensidade de espalhamento dos modelos das diferentes abas para exibir no gráfico da interface'''
        total = []
        for j in range(len(self.ff_file)):
            total.append(np.zeros(len(self.ff_file[j].q)))
            for i in range(len(self.intensidade)):
                if len(self.intensidade[i])==len(total[j]) and self.relative_file[i].currentIndex()-1 == j:
                    if i==0:
                        total[j] = np.zeros(len(self.intensidade[i]))
                    #print('j, {}, i: {}, total[j]: {}'.format(j, i, total[j]))
                    total[j] = np.add(total[j], self.intensidade[i]*10**self.lmf_params['scale'+'_comma_'+'Model'+'_comma_'+str(i)].value)
        return total
    
    #########################   FORM FACTOR FUNCTIONS    ###########################    
    def prepare_ff_widgets(self):
        '''Essa função cria os arrays onde serão salvos os dicionários de cada widget necessário nas diferentes abas'''
        self.intensidade = []
        self.expr = []
        #self.lmf_params = []
        self.ff_param_groupbox = []
        self.ff_layout = []
        self.ffbox = []
        self.addtab_button = []
        self.removetab_button = []
        self.methodbox = []
        self.sf_Ref_textbox = []
        self.scale = []
        self.relative_file = []
    
    def addFormFactors(self): 
        '''Esta função adiciona a opção de simular fatores de forma'''

        self.q = np.logspace(-3.2,0, np.random.randint(500, high=1000)) #vetor q
        print(len(self.q))

        self.prepare_ff_widgets()        

        self.add_ff_tab_function()
        #calcula o primeiro fator de forma para plotar no grafico
        #self.intensidade[0] = self.ffs[0].intensidade
        self.intensidade_plot = self.sum_intensidade()
 
    def remove_ff_widgets(self):
        '''Essa função apaga os widgets quando uma aba é removida'''
        index = self.ff_tabs.currentIndex()

        self.relative_file.pop(index)
        self.ff_widgets_added.pop(index)
        self.dist_widgets_added.pop(index)
        self.sf_widgets_added.pop(index)

        self.tabs.pop(index) 

        self.ff_param_groupbox.pop(index)
        self.ff_layout.pop(index)

        self.intensidade.pop(index)
        self.method.pop(index)

        self.ffbox.pop(index)

        self.addtab_button.pop(index)
        self.removetab_button.pop(index)
        self.methodbox.pop(index)

        self.gb_layout.pop(index)

        self.ffs.pop(index)
        self.remove_ff_pars_widgets()
    
    def remove_ff_pars_widgets(self):
        '''Essa função apaga os widgets dos sliders quando uma aba é removida'''
        index = self.ff_tabs.currentIndex()

        self.slider['Model'].pop(index)
        self.param_text['Model'].pop(index)
        self.text_edit_min['Model'].pop(index)
        self.text_edit_max['Model'].pop(index)
        self.text_edit_value['Model'].pop(index)
        self.param_vary_cb['Model'].pop(index)
        self.param_expr_textbox['Model'].pop(index)
    
    def add_ff_tab_function(self):
        '''Essa função adiciona uma nova aba de fatores de forma'''
        index = self.ff_tabs.currentIndex()

        if len(self.tabs)>0:
            fforma = self.ffbox[index].currentIndex()
        else:
            fforma = 1

        self.ff_widgets_added.append(False)
        self.dist_widgets_added.append(False)
        self.sf_widgets_added.append(False)

        self.tabs.append(QtWidgets.QWidget()) 
        self.tabs[-1].layout = QtWidgets.QVBoxLayout()  

        self.ff_tabs.addTab(self.tabs[-1], 'Tab {}'.format(len(self.tabs)))   
        self.ff_tabs.setCurrentIndex(len(self.tabs)-1)
        index = self.ff_tabs.currentIndex()

        self.ff_param_groupbox.append(QtWidgets.QGroupBox('Form Factor'))
        self.ff_layout.append(QtWidgets.QGridLayout())

        self.method.append(0)

        self.addNeededWidgets()

        ff = self.which_ff(fforma)
        ff.setq(self.q)
        ff.ff()
        self.ffs.append(ff)

        self.intensidade.append(ff.intensidade)

        self.ffbox.append(QtWidgets.QComboBox()) #caixa com as opções para mudar o fator de forma
        self.ffbox[-1].addItems(f.all)
        self.ffbox[-1].setCurrentIndex(fforma)
        self.ffbox[-1].currentIndexChanged.connect(self.ffChanged)

        self.addtab_button.append(QtWidgets.QPushButton('+'))
        self.addtab_button[-1].setMaximumWidth(40)
        self.removetab_button.append(QtWidgets.QPushButton('-'))
        self.removetab_button[-1].setMaximumWidth(40)
        if index == 0:
            self.removetab_button[-1].setEnabled(False)

        self.addtab_button[-1].clicked.connect(self.add_ff_tab_function)
        self.removetab_button[-1].clicked.connect(self.removetab_function)

        self.methodbox.append(QtWidgets.QComboBox()) #caixa com as opções para mudar o método dos inputs dos fatores de forma
        self.methodbox[-1].currentIndexChanged.connect(self.ff_methodChanged)   

        self.relative_file.append(QtWidgets.QComboBox())  
        self.relative_file[-1].setToolTip('Selecione o Arquivo')
        self.add_relative_file_cb_options()

        self.position['Model'][-1]=5 #formatação visual
        self.addparam('scale', -10, 10, 0, 'Scale', self.ff_layout[-1]) #adiciona o slider de escala

        #variavel 'self.intensidade_plot' é a curva que será plotada. 
        #É o fator de forma calculado é multiplicado por uma escala e somado a um background

        self.ff_addparam() #adiciona os sliders dos parâmetros necessários para calcular o fator de forma do objeto escolhido
        self.ff_layout[-1].addWidget(self.ffbox[-1], 0,5)
        self.ff_layout[-1].addWidget(self.methodbox[-1], 1,5) 
        self.ff_layout[-1].addWidget(self.relative_file[-1], 2, 5)
        self.ff_layout[-1].addWidget(self.addtab_button[-1], 0, 8)
        self.ff_layout[-1].addWidget(self.removetab_button[-1], 0, 9)
        #self.adj_graph() #ajusta os limites do grafico
        self.position['Model'][-1]+=1 #formatação visual

        self.ff_param_groupbox[-1].setLayout(self.ff_layout[-1])
        self.ff_param_groupbox[-1].setMaximumHeight(290)
        self.ff_param_groupbox[-1].setMinimumHeight(290)

        self.gb_layout.append(QtWidgets.QVBoxLayout())

        #self.sf_Ref_textbox.append('')
        self.expr.append('')

        self.gb_layout[-1].addWidget(self.ff_param_groupbox[-1])

        self.tabs[-1].setLayout(self.gb_layout[-1])  

        self.scale.append(10**self.slider['Model'][-1]['scale'].value())

        self.relative_file[-1].setCurrentIndex(0)

        self.relative_file[-1].currentIndexChanged.connect(lambda: self.relative_file_changed(index))

        self.ff_widgets_added[-1] = True #variável de controle para saber se os widgets daquela aba ja foram adiconados.
        
        #TODO essa parte eu não sei se vale a pena mudar ou manter assim
        if len(self.tabs)>1:
            if self.dist_added:
                self.add_dist_tab()
                self.dist_widgets_added[index] = True
            if self.sf_added:
                self.add_sf_tab()
                self.sf_widgets_added[index] = True
            self.create_lmf_params()

        self.valuechange('scale')

        return self.gb_layout[-1]

    def remove_ff_pars(self):
        '''Esta função remove apenas os sliders referentes aos parâmetros necessários para calcular o fator de forma do objeto escolhido'''
        index = self.ff_tabs.currentIndex()

        for i in (self.ffs[index].func_pars):
            self.slider['Model'][index][i].deleteLater()
            self.param_text['Model'][index][i].deleteLater()
            self.param_vary_cb['Model'][index][i].deleteLater()
            self.text_edit_min['Model'][index][i].deleteLater()
            self.text_edit_max['Model'][index][i].deleteLater()
            self.text_edit_value['Model'][index][i].deleteLater()
            self.param_expr_textbox['Model'][index][i].deleteLater()
            self.slider['Model'][index].pop(i)

    def which_ff(self, fforma=0):
        '''Essa função retorna o objeto do modelo escolhido.
        A variável "fforma" informa qual o modelo escolhido.'''
        ff = f.dic[fforma]
        func = ff()
        return func
       
    def ffChanged(self): 
        '''Esta função muda os sliders. É chamada ao trocar o formato do objeto cujo fator de forma esta sendo simulado'''
        index = self.ff_tabs.currentIndex()

        #Trocas as variáveis de controle que informam se os widgets ja foram adiconados para 'False', pois eles serão trocados
        self.ff_widgets_added[index] = False
        if self.dist_added:
            self.dist_widgets_added[index] = False
        if self.sf_added:
            self.sf_widgets_added[index] = False

        #Desativa as checkboxs das opções de aplicar distribuição e fator de estrutura
        if self.dist_added:
            self.dist_cb[index].setChecked(False)
        if self.sf_added:
            self.sf_cb[index].setChecked(False)

        q = self.ffs[index].q

        if bool(self.slider['Model'][index]): #verifica se havia sliders do antigo objeto, se sim, os remove
            self.remove_ff_pars()
        fforma = self.ffbox[index].currentIndex() #verifica qual o formato do objeto
        ff = self.which_ff(fforma) 
        self.ffs[index] = ff #salva o novo objeto fator de forma
        if self.relative_file[index].currentIndex()==0: #se o fator de formar não for relacionado à nenhum arquivo de dados, para o cálculo e usado um 'q' genérico.
            self.ffs[index].setq(self.q)
        else: #se não é usado o 'q' do arquivo de dados
            self.ffs[index].setq(q)

        if self.dist_added:
            if self.ffs[index].name == 'Background':
                self.dist_cb[index].setEnabled(False)
            else:
                self.dist_cb[index].setEnabled(True)
        if self.sf_added:
            if self.ffs[index].name == 'Background':
                self.sf_cb[index].setEnabled(False)
            else:
                self.sf_cb[index].setEnabled(True)

        self.methodbox[index].clear() #limpa as opções da caixa 

        #self.position=6 #formatação visual
        self.ff_addparam() #adiciona os sliders dos parâmetros
        
        if self.methodbox[index].currentIndex() >= 0:
            self.param_box[index].clear()
            self.param_box[index].addItems(self.ffs[index].dist_params[self.methodbox[index].currentIndex()])

        self.ff_widgets_added[index] = True
        if self.dist_added:
            self.dist_widgets_added[index] = True
        if self.sf_added:
            self.sf_widgets_added[index] = True
        self.create_lmf_params() #cria os parâmetros do lmfit

        self.valuechange('scale')
        
    def ff_addparam(self): 
        '''Esta função adiciona os sliders dos parâmetros necessário para o cálculo do fator de forma do objeto escolhido'''
        index = self.ff_tabs.currentIndex()

        for i in self.ffs[index].func_pars: #cria um slider para cada parâmetros do modelo definido dentro do dicionário 'func_pars'
            self.ff_layout[index] = self.addparam(i, self.ffs[index].func_pars[i]['min'], self.ffs[index].func_pars[i]['max'], 
                                                                    self.ffs[index].func_pars[i]['value'], self.ffs[index].func_pars[i]['symbol'], self.ff_layout[index])
        self.methodbox[index].addItems(self.ffs[index].methods) #adicona as opções de variáveis de input na combobox da interface

        self.ff_param_groupbox[index].setLayout(self.ff_layout[index]) #formatação visual
        #self.right_layout.addWidget(self.ff_param_groupbox, 0)
    
    def ff_methodChanged(self): 
        '''Esta variável muda o método do input dos parâmetros. '''
        index = self.ff_tabs.currentIndex()

        method = self.methodbox[index].currentIndex()
        if method>=0:# and self.ff_widgets_added[index]:
            if len(self.ffs[index].changeparam)>1: #verifica se há mais de um metodo
                for i in range(len(self.ffs[index].methods)):
                    if method == i: #verifica quais a variáveis de input escolhidas
                        self.method[index] = i 
                        self.ffs[index].value_methods(i, self.slider['Model'][index]) #converte os valores de input entre o métodos
                self.ff_changeParamFunc(self.ffs[index].changeparam[0], self.ffs[index].func) #muda a função aplicada ao valor do slider
                self.ff_changeparam(self.ffs[index].changeparam[0], self.ffs[index].changeparam[1],
                                    self.ffs[index].value['min'], self.ffs[index].value['max'], 
                                    self.ffs[index].value['value']) #troca o parametro
                                
                if self.dist_added and self.dist_widgets_added[index]:
                    #troca as opções de variáveis vísiveis na caixa de selceção da distribuição para as novas variáveis
                    self.param_box[index].clear()
                    self.param_box[index].addItems(self.ffs[index].dist_params[self.methodbox[index].currentIndex()])
        #self.valuechange('scale')   

    def ff_changeParamFunc(self, name, eq):
        '''Esta função muda a função aplicada ao valor do slider
        A variável "name" é uma string informando qual parâmetro do modelo mudou de equação.
        A variável "eq" é uma string com a equação fornecida.'''
        self.slider['Model'][self.ff_tabs.currentIndex()][name].func_value = eq
        self.valuechange(name)

    def ff_removeAllPar(self):
        '''Esta função remove todos os sliders do layout'''
        index = self.ff_tabs.currentIndex()

        for i in (self.slider[index]): #remove os sliders
            self.slider['Model'][index][i].deleteLater()
            self.param_text['Model'][index][i].deleteLater()
            self.param_vary_cb['Model'][index][i].deleteLater()
            self.text_edit_min['Model'][index][i].deleteLater()
            self.text_edit_max['Model'][index][i].deleteLater()
            self.text_edit_value['Model'][index][i].deleteLater()
            self.param_expr_textbox['Model'][index][i].deleteLater()
            self.slider['Model'][index].clear()
        #apaga os dicionários
        self.slider['Model'][index] = {}
        self.param_text['Model'][index] = {}
        self.param_vary_cb['Model'][index] = {}
        self.text_edit_min['Model'][index] = {}
        self.text_edit_max['Model'][index] = {}
        self.text_edit_value['Model'][index] = {}
        self.param_expr_textbox['Model'][index] = {}

    def ff_changeparam(self, name, newname, newmin, newmax, newvalue): 
        '''Esta função muda o parâmetro do slider e define os novos valores máximos, mínimos e atuais do slider. 
        Ela é usada quando o método de input do parâmetro é trocado na interface.
        A variável "name" é uma string informando qual parâmetro do modelo mudou de valor.'''
        index = self.ff_tabs.currentIndex()

        self.param_text['Model'][index][name].setText(str(newname))
        self.text_edit_value['Model'][index][name].setText(str(newvalue))
        self.text_edit_min['Model'][index][name].setText(str(round(newmin, 2)))
        self.text_edit_max['Model'][index][name].setText(str(round(newmax, 2)))

        self.slider['Model'][index][name].setMinimum(round(newmin, 2))
        self.slider['Model'][index][name].setMaximum(round(newmax, 2))
        self.slider['Model'][index][name].setValue(newvalue)
        self.valuechange(name)

    #######################   FF TABS FUNCTION    ##########################
    def removetab_function(self):
        '''Essa função remove uma aba e apaga seus widgets'''
        tab = self.ff_tabs.currentIndex()
        self.remove_ff_widgets()
        if self.dist_added:
            self.remove_dist_widgets()
        if self.sf_added:
            self.remove_sf_widgets()
        self.ff_tabs.removeTab(self.ff_tabs.currentIndex())
        for i in range(tab, len(self.tabs)):
            self.ff_tabs.setTabText(i, 'Tab {}'.format(i+1))
        if tab < len(self.tabs):
            self.ff_tabs.setCurrentIndex(tab)
        else:
            self.ff_tabs.setCurrentIndex(len(self.tabs)-1)
        self.create_lmf_params()
        self.valuechange('scale')
    
    ########################   RAIO DE GIRO    #################################
    def addRg(self): 
        '''Esta função adiciona a opção de calcular o raio de giro por guinier de um dado importado'''
        self.rg_layout = QtWidgets.QGridLayout()
        self.rg_groupbox = QtWidgets.QGroupBox('Raio de Giro')
        self.raiodegiro = f.RaioDeGiro(self.ffs[self.ff_tabs.currentIndex()]) #cria o objeto raio de giro
        self.qRg_text = QtWidgets.QLabel() 
        self.qRg_text.setText('qRg: ')
        self.qRg_text.setStyleSheet(self.SS())
        self.Rg_text = QtWidgets.QLabel()
        self.Rg_text.setStyleSheet(self.SS())
        self.Rg_text.setText('Rg: ')
        self.qRg_textbox = QtWidgets.QLineEdit()
        self.qRg_textbox.setMaximumWidth(80)
        self.qRg_textbox.setText('1')
        self.qRg_textbox.returnPressed.connect(self.raioDeGiro)
        self.Rg_textbox = QtWidgets.QLineEdit()
        self.Rg_textbox.setMaximumWidth(80)
        self.Rg_textbox.setEnabled(False)
        #self.Rg_textbox.setText(str(round(self.ffs[self.ff_tabs.currentIndex()].Rg, 2)))
        self.rg_layout.addWidget(self.qRg_text, 0, 0)
        self.rg_layout.addWidget(self.Rg_text, 1, 0)
        self.rg_layout.addWidget(self.qRg_textbox, 0, 1)
        self.rg_layout.addWidget(self.Rg_textbox, 1, 1)
        self.rg_layout.setAlignment(QtCore.Qt.AlignLeft)
        self.rg_groupbox.setLayout(self.rg_layout)
        self.rg_groupbox.setMaximumHeight(90)
        self.left_layout.addWidget(self.rg_groupbox, 4)

    def raioDeGiro(self):
        '''Esta função calcula o raio de giro do dado importado e exibe o valor na caixa de texto'''
        qRg = float(self.qRg_textbox.text()) #verifica o valor máximo de q*Rg
        self.raiodegiro.set_qRg(qRg) 
        if self.ff_file[-1]: #verifica de algum dado foi importado
            self.raiodegiro.setff(self.ff_file[-1]) #calcula o raio de giro
            if self.ff_file[-1].Rg:
                self.Rg_textbox.setText(str(round(self.ff_file[-1].Rg, 2))) #exibe o valor calculado

    ################      FIT FUNCTIONS       #############################
    def fit_apply_button_clicked(self):
        #self.create_lmf_params()
        self.x_fit = [] #array com os valores de x dos dados para o fit
        self.y_fit = [] #array com os valores de y dos dados para o fit
        self.y_error_fit = [] #array com os valores de erro de y dos dados para o fit
        #self.index_fit = []
        #print('len: {}'.format(len(self.lmf_params_tabs)))
        self.create_lmf_params()

        for i in range(len(self.ff_file)): #for para passar por todos os dados a serem fitados
            if self.file_fit_cb[i].isChecked(): #verifica se o dado está selecionado na interface gráfica para ser fitado
                self.x_fit.append(self.ff_file[i].q[self.indmin_f[i]:self.indmax_f[i]]) #adiciona o array com os valores de x do dado na variável x_fit
                self.y_fit.append(self.ff_file[i].intensidade[self.indmin_f[i]:self.indmax_f[i]]) #mesma coisa para os valore de y
                if self.ff_file[i].yerr is None:
                    self.y_error_fit.append(self.apply_eps(self.ff_file[i])[self.indmin_f[i]:self.indmax_f[i]]) #mesma coisa para os valores do erro de y
                else:
                    self.y_error_fit.append(self.ff_file[i].yerr[self.indmin_f[i]:self.indmax_f[i]])
                #self.index_fit.append(i)
            else:
                pass

        self.out = (lmf.minimize(self.res, self.lmf_params, args=(self.x_fit,), kws={'data':self.y_fit, 'eps':self.y_error_fit}))#, 'index':self.index_fit}))
        self.fit_report_popup()
        self.ax.legend(loc='upper left', prop={'size':8})

    def fit_upd_report(self):
        self.fit_report_text_widget[-1].setText(lmf.fit_report(self.out).replace('_comma_', ''))
        #print((lmf.fit_report(self.out)))
    
    def res(self, pars, x, data=None, eps=None):#, index=None):
        ndata = len(x) #número de curvas a serem fitadas
        model = [] 
            
        last = [] #array com o número de onde cada parâmetro é 
        scale = [] #array com o valor da escala para cada parâmetro
        parvals = pars.valuesdict() #dicionário com os valores de cada parâmetro
        parvals_keys = list(parvals.keys()) #array com as chaves do diconário "parvals"
        for i in parvals: #percorrendo cada item do dicionário
            #print(i)
            exec(i+'= parvals[i]') #cria uma variável para cada parâmetro cujo nome é a chave do dicionário
            i_ = i.split('_comma_')
            last.append(i_[2])
            #print(i, i[:5])
            if i_[0] == 'scale': #verifica se o parâmetro é referente à escala
                #print('entrei, i={}'.format(i))
                scale.append(parvals[i]) #adiciona o valor no array de escala
         
        for i in range(ndata): #percorre por todas as curvas a serem fitadas
            funcpars = [] #array com os dricionário "funcpars" com os parâmetros de cada fator de forma 
            funcpars_sf = [] #o mesmo para o fator de estrutura
            model.append(np.zeros(len(x[i]))) #adiciona uma coluna de zeros com o tamanho do número de pontos dos dados
            for j in range(len(self.relative_file)): #percorre por cada aba
                if i == self.relative_file[j].currentIndex()-1: #verifica se a aba esssa associada com o dado atual
                    funcpars.append(self.ffs[j].func_pars) #adiciona o funcionário "funcpars"
                    #if self.sf_cb[j].isChecked():
                    funcpars_sf.append(self.sf[j].func_pars) #adiciona o dicionário "funcpars" do fator de estrutura
                    for l in range(len(parvals_keys)): #percorre por cada parâmetro
                            if float(last[l]) == j: #verifica se o parâmetro é relativo a tabela atual
                                index_ = parvals_keys[l].split('_comma_')[0]
                                #index_.pop(-1)
                                #index_ = ''.join(index_) #the index_ variable is the parameter name without the number relative to the tab
                                #print(index_)
                                #print('funcpars[j]: ', funcpars[j])
                                if index_ != 'scale' and index_ != 'dist': #verifica se o parâmetro é a variável "scale" ou "dist", pois elas não fazer parte dos dicionários de parâmetro nem do fator de estrutura nem do fator de forma
                                    if index_ in list(self.ffs[j].func_pars.keys()): #verifica se o parâmetro é da função do fator de forma
                                        funcpars[-1][index_]['value'] = parvals[parvals_keys[l]] #adiciona o valor ao dicionário "funcpars"
                                    elif index_ in list(self.sf[j].func_pars.keys()): #verifica se o parâmetro é da função do fator de estrutua
                                        funcpars_sf[-1][index_]['value'] = parvals[parvals_keys[l]] #adiciona o valor ao dicionário "funcpars_sf"
                    #print('funcpars: ', funcpars)
                    #print('len scale: {}, len ffs: {}, len funcpars: {}'.format(len(scale), len(self.ffs), len(funcpars)))

                    #verifica se é necessário aplicar a distribuição de tamanhos ou o fator de estrutra e chama a função apropriada
                    if self.dist_cb[j].isChecked() and not(self.sf_cb[j].isChecked()):
                        #print(funcpars)
                        model[i]+=self.apply_dist(j, x[i], parvals['dist'+'_comma_'+'Model'+'_comma_'+str(j)], funcpars[-1], False)*(10**scale[j])
                    elif not(self.dist_cb[j].isChecked()) and self.sf_cb[j].isChecked():
                        model[i]+=self.apply_sf(None, j, x[i], funcpars[-1], funcpars_sf[-1], False)*(10**scale[j])
                    elif self.dist_cb[j].isChecked() and self.sf_cb[j].isChecked():
                        model[i]+=self.apply_sf_dist(None, j, x[i], parvals['dist'+'_comma_'+'Model'+'_comma_'+str(j)], funcpars[-1], funcpars_sf[-1], False)*(10**scale[j])
                    else:
                        #print((self.ffs[j].ff3(x, funcpars[-1])*(10**scale[-1])).shape)
                        model[i]+=(self.ffs[j].ff3(x[i], funcpars[-1])*(10**scale[j]))
        if data is None:
            return model
        if eps is None:
            return (self.flatten(model) - self.flatten(data))
        return ((self.flatten(model)-self.flatten(data)) / self.flatten(eps))
    
    def apply_eps(self, data):
        delta_error = np.logspace(0.001, 2, len(data.q))
        return abs(data.intensidade*delta_error)
    
    def fit_report_popup(self):
        '''Essa função cria um popup sobre o sucesso do fit'''
        self.fit_msg = QtWidgets.QMessageBox() 
        self.fit_msg.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.fit_msg.setIcon(QtWidgets.QMessageBox.Question)
        self.fit_msg.setWindowTitle('Fit Report')
        success = self.out.success
        if success: #mensagem de sucesso no fit
            self.fit_msg.setText('Success! \n Accept Fit? \n Click on details to see the complete fit report.')
            #message = lmf.fit_report(self.out)
            #self.fit_msg.setDetailedText(message)
        elif success == False: #mensagem de erro no fit
            self.fit_msg.setText("Couldn't to perform the fit.") #mensagem
            for i in range(len(self.relative_file)): #percorre pela abas de fatores de forma
                self.ffs[i].func_pars = self.old_funcpars_ff[i] #garante que os valores dos parâmetros nos objetos da classe "Form Factors" são os mesmo que antes da tentativa de fit
                if self.sf_added and self.sf_cb[i].isChecked(): #verifica se as apções de fator de estrutura foram adicionadas e se estão habilitadas na aba atual
                    self.sf[i].func_pars = self.old_funcpars_sf[i] #garante que os valores dos parâmetros dos objetos da classe "Structure Factor"são os mesmos que antes da tentativa de fit
        elif success == None: #mensagem caso nenhum dado tenha sido selecionado para o fit
            self.fit_msg.setText("No data selected to fit.")
                             
        self.fit_msg.setStandardButtons(QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.Ok)
        self.fit_msg.buttonClicked.connect(self.accept_fit)
        self.fit_msg.exec_()

    def accept_fit(self, button):
        '''Função executada quando qualquer botão do popup sobre o fit é clicado'''
        button = self.fit_msg.clickedButton() #botão clicado
        sb = self.fit_msg.standardButton(button)
        if sb == QtWidgets.QMessageBox.Ok: #verifica se o botão clicado foi o botão "Ok"
            self.add_fit_tab() #cria uma nova aba de fit report
            model = self.res(self.out.params, self.x_fit) #calcula a curva fitada e salva na variável "Model"
            index_model = 0
            for i in range(len(self.ff_file)): #percorre por cada arquivo de dado
                #! talvez o indice i esteja errado, verificar
                if self.file_fit_cb[i].isChecked(): #verifica se o dado foi selecionado para ser fitado
                    self._line3[i].set_data(self.x_fit[i], model[index_model])
                    self._line3[i].set_visible(True)
                    index_model += 1
            self._line3[-1].figure.canvas.draw()  #cria e desenha a curva fitada no gráfico da interface
            self.fit_upd_report()
            #self.fit_save_data_button.setEnabled(True)
            #self.fit_plot_res_button.setEnabled(True)
            
        elif sb == QtWidgets.QMessageBox.Cancel: #verifica se o botão clicado foi o botão "Cancel"
            for i in range(len(self.relative_file)): #percorre por cada aba de fator de forma
                #garante que os valores dos parâmetro nos objetos das classes de fator de forma e de estrutura são os mesmos que antes da tentativa de fit
                self.ffs[i].func_pars = self.old_funcpars_ff[i] 
                self.sf[i].func_pars = self.old_funcpars_sf[i]

    def create_lmf_params(self):
        '''Essa função cria os obejtos da classe "Parameters" da biblioteca lmFit'''
        self.lmf_params = NewLmfParameters(self) #cria o objeto
        self.lmf_params._asteval.symtable["np"] = np
        #print(type(self.lmf_params))
        self.old_funcpars_ff = [] #array onde serão salvos os valores de cada dicionário "funcpars" de antes do fit para cada aba de fator de forma
        self.old_funcpars_sf = [] #mesma coisa para os valores do fator de estrutua
        _expr = []

        for tab_ in range(self.adj_tabs.count()):
            tab = self.adj_tabs.tabText(tab_)  
            if tab == 'Model':  
                for j in range(len(self.relative_file)): #percorre por cada aba de fator de forma
                    self.old_funcpars_ff.append(self.ffs[j].func_pars)
                    if self.sf_widgets_added[j] and self.sf_cb[j].isChecked():
                        self.old_funcpars_sf.append(self.sf[j].func_pars)
                        #if i == self.relative_file[j].currentIndex()-1:
                    for k in self.slider[tab][j]:
                        #print('_expr_: ', _expr_)
                        _expr_ = self.param_expr_textbox[tab][j][k].text()
                        if _expr_ == '':
                            _expr.append(None)
                        else:
                            _expr.append(_expr_)
                        if self.slider[tab][j][k].isEnabled(): #verifica se o slider está ativo
                            #print(i)
                            #print('min: {}, max: {}'.format(self.slider[j][k].minimum(), self.slider[j][k].maximum()))
                            
                            #se estiver cria o parâmetro com o valor do slider
                            #a string "_comma_" está sendo usada para separa a parte do nome do parâmetro que é referente ao nome do parâmetro, à aba geral da interface e à aba do fator de forma
                            self.lmf_params.add(k+'_comma_'+tab+'_comma_'+str(j), value=self.slider[tab][j][k].value(), vary=self.param_vary_cb[tab][j][k].isChecked(), min=self.slider[tab][j][k].minimum(),
                                                max=self.slider[tab][j][k].maximum())
                        else:
                            #print(k)
                            self.lmf_params.add(k+'_comma_'+tab+'_comma_'+str(j), value=None, vary=False, min=self.slider[tab][j][k].minimum(), max=self.slider[tab][j][k].maximum())
            elif tab == 'Peaks':
                for j in range(len(self.peaks_tabs_array)): #percorre por cada aba de dos picos
                    #print('Peaks tabs: ', j)
                    for k in self.slider[tab][j]: #percorre cada parâmetro de cada aba
                        _expr_ = self.param_expr_textbox[tab][j][k].text()

                        #adicona a expressão ao parâmetro
                        if _expr_ == '':
                            _expr.append(None)
                        else:
                            _expr.append(_expr_)
                        
                        if self.slider[tab][j][k].isEnabled(): #verifica se o slider está ativo
                            #print(i)
                            #print('min: {}, max: {}'.format(self.slider[j][k].minimum(), self.slider[j][k].maximum()))
                            
                            #se estiver cria o parâmetro com o valor do slider
                            self.lmf_params.add(k+'_comma_'+tab+'_comma_'+str(j), value=self.slider[tab][j][k].value(), vary=self.param_vary_cb[tab][j][k].isChecked(), min=self.slider[tab][j][k].minimum(),
                                                max=self.slider[tab][j][k].maximum())
                        else:
                            #se não cria o parâmetro com valor nulo
                            self.lmf_params.add(k+'_comma_'+tab+'_comma_'+str(j), value=None, vary=False, min=self.slider[tab][j][k].minimum(), max=self.slider[tab][j][k].maximum())
            
            for i in self.lmf_params.valuesdict():
                #remove os arrays "_comma_", "Model" e "Peaks" do nome dos parâmetros e salva na variável "name", deixando apenas o nome do parâmetro ('raio' por exemplo e o número da aba em que ele está)
                name = i.replace('_comma_', '')
                name = name.replace('Model', '')
                name = name.replace('Peaks', '')
                self.lmf_params._asteval.symtable[name] = self.lmf_params[i].value #adiciona esse "novo" novo na symbol table do lmfit, para que esse parâmetro possa ser usado de forma mais simples na caixa de expressões da interface gráfica
            self.lmf_add_expr(_expr) #adiciona a expressão matemática a todos os parâmetros

    def lmf_add_expr(self, _expr):
        '''Essa função adciona as expressões matemáticas a todos os parâmetros do lmfit.
        A variável "_expr" é um array com as expreções de cada parâmetro'''
        #! Vendo agora acredito que a maneira que está feita aqui não seja a melhor. Pode ser que eu esteja atribuindo expressões à parâmetros errados. Preciso verificar.
        #! Se for preciso mudar essa função também será preciso alterar a função "create_lmf_params".
        counter = 0
        for i in self.lmf_params.keys(): #percorre todos os parâmetros
            self.lmf_params[i].expr = _expr[counter] #adiociona a expressão
            counter+=1
                
    def add_fit_tab(self):
        '''Essa função adiciona uma aba com os dados do fit realizado'''
        self.all_fit_record.append(self.out)
        self.fit_record.append(self.out)

        self.fit_layout.append(QtWidgets.QGridLayout())
        #self.fit_groupbox.append(QtWidgets.QGroupBox())
        self.fit_tabs_array.append(QtWidgets.QWidget())
        self.fit_save_data_button.append(QtWidgets.QPushButton())
        self.fit_plot_res_button.append(QtWidgets.QPushButton('Plot Residual'))
        self.fit_remove_tab_button.append(QtWidgets.QPushButton())
        self.fit_send_params_button.append(QtWidgets.QPushButton())
        #self.fit_tabs_array[-1].layout = QtWidgets.QGridLayout()
        self.fit_report_text_widget.append(QtWidgets.QTextEdit('No fit performed for this file.'))
        self.fit_report_text_widget[-1].setReadOnly(True)
        self.fit_report_text_widget[-1].setMaximumHeight(400)
        self.fit_save_data_button[-1].setIcon(qta.icon('fa.save'))
        self.fit_save_data_button[-1].setToolTip('Save fit data')
        self.fit_remove_tab_button[-1].setIcon(qta.icon('ri.delete-bin-5-fill'))
        self.fit_remove_tab_button[-1].setToolTip('Remove this fit')
        self.fit_send_params_button[-1].setIcon(qta.icon('fa.check-circle-o'))
        self.fit_send_params_button[-1].setToolTip('Send fit data to sliders')
        self.fit_save_data_button[-1].clicked.connect(self.fit_save_data_button_clicked)
        self.fit_plot_res_button[-1].clicked.connect(self.add_fit_res_graph)
        self.fit_remove_tab_button[-1].clicked.connect(self.remove_fit_tab)
        self.fit_send_params_button[-1].clicked.connect(self.fit_send_params)
        self.fit_layout[-1].addWidget(self.fit_save_data_button[-1], 0, 0)
        self.fit_layout[-1].addWidget(self.fit_plot_res_button[-1], 0, 1)
        self.fit_layout[-1].addWidget(self.fit_send_params_button[-1], 0, 2)
        self.fit_layout[-1].addWidget(self.fit_remove_tab_button[-1], 0, 3)
        self.fit_layout[-1].addWidget(self.fit_report_text_widget[-1], 1, 0, 1, 4)
        self.fit_tabs_array[-1].setLayout(self.fit_layout[-1])
        self.fit_tabs.addTab(self.fit_tabs_array[-1], 'Fit '+str(len(self.fit_record)))
        self.fit_final_layout.addWidget(self.fit_tabs)   
        self.fit_tabs.setCurrentIndex(len(self.fit_tabs_array)-1)         

    def fit_send_params(self):
        '''Essa função envia os valores dos parâmetro encontrados no fit para os sliders da interface.'''
        #! Ela pode ter algum problema, mas não consigo encontrar. Quando um fit que envolve mais de uma aba é realizado e essa função é utilizada ele envia os parâmetros corretamente para apenas uma aba.
        #! Nos parâmetros das outras abas a posição do slider é alterada mas seu valor não. Já alterei essa fução várias vezes a ainda não consegui encontrar a causa do problema.
        for i in range(len(self.ff_file)): #percorre todos os arquivos de dados
            if self.file_fit_cb[i].isChecked(): #verifica se estão selecionados para serem fitados
                self._line3[i].set_visible(False)
                for tab_ in range(self.adj_tabs.count()): #percorres as abas principais da iterface ("Data", "Model" e "Peaks")
                    tab = self.adj_tabs.tabText(tab_) #salva o nome da aba na variável tab
                    if tab == 'Model': #verifica se a aba é a "Model" (aba dos fatores de forma)
                        for j in range(len(self.relative_file)): #percorre todas a abas de fatores de forma (são as abas dentro da aba "Model")
                            if i == self.relative_file[j].currentIndex()-1: #verifica se o dado relativo a essa aba é o dado atual
                                for k in self.slider[tab][j]: #percorre sliders dos parâmetros da aba
                                    if self.slider[tab][j][k].isEnabled(): #verifica se o slider do parâmetro está ativo
                                        self.slider[tab][j][k].setValue(self.fit_record[self.fit_tabs.currentIndex()].params.valuesdict()[k+'_comma_'+tab+'_comma_'+str(j)]) #altera o valor do slider para o valor encontrado no fit
                                        #self.slider[tab][j][k].sliderPosition = int(self.fit_record[self.fit_tabs.currentIndex()].params.valuesdict()[k+'_comma_'+tab+'_comma_'+str(j)])
                                        #self.slider[tab][j][k].update()
                                        #self.slider[tab][j][k].repaint()
                                        self.valuechange(k, j, False, 'Model', False) #chama a função valuechange para atualizar o valor dos parâmetros que dependam via as expressões do que acabou de ser alterado
                                self.valuechange(k, j, True, 'Model', True) #recalcula a intensidade de espalhamento e exibe a nova curva no gráfico da interface
        
    def prepare_fit_widgets(self):
        '''Essa função prepara os widgets das abas onde são exibidos os dados dos ajustes realizados'''
        #self.fit_groupbox = []
        self.fit_layout = []
        self.fit_tabs = QtWidgets.QTabWidget()
        self.fit_tabs_array = []
        self.fit_report_text_widget = []
        self.fit_save_data_button = []
        self.fit_plot_res_button = []
        self.fit_remove_tab_button = []
        self.fit_send_params_button = []
        self.fit_record = []

        self.fit_final_layout = QtWidgets.QVBoxLayout()        
        self.fit_apply_button = QtWidgets.QPushButton('Fit Data')
        self.fit_apply_button.clicked.connect(self.fit_apply_button_clicked)
        self.fit_final_layout.addWidget(self.fit_apply_button)
        self.fit_final_layout.setAlignment(QtCore.Qt.AlignTop)
        self.adj_tabs_array[2].setLayout(self.fit_final_layout)

    def remove_fit_tab(self):
        '''Essa função remove uma aba de dados do fit'''
        index = self.fit_tabs.currentIndex()
        self.remove_fit_widgets(index)
        self.fit_tabs.removeTab(index)
        #self._line3.pop(index)

    def remove_fit_widgets(self, index):
        '''Essa função remove os widgets das abas onde são exibidos os dados dos ajustes realizados'''
        self.fit_layout.pop(index)
        self.fit_report_text_widget.pop(index)
        self.fit_tabs_array.pop(index)
        self.fit_save_data_button.pop(index)
        self.fit_plot_res_button.pop(index)
        self.fit_remove_tab_button.pop(index)
        self.fit_send_params_button.pop(index)
        self.fit_record.pop(index)
    
    def add_fit_res_graph(self):
        '''Essa função exibie o gráfico de resíduo do fit.'''
        self.fit_sim_figure, self.fit_ax = plt.subplots()
        initial = 0
        end = 0
        for i in range(len(self.ff_file)):
            end += len(self.ff_file[i].q)
            if self.file_fit_cb[i].isChecked():
                self.fit_ax.plot(self.ff_file[i].q, self.out.residual[initial:end], label=self.label_s[i])
            initial = end
        self.fit_ax.set_xscale('log') #muda a escala do eixo x do gráfico
        self.fit_ax.set_yscale('log') #muda a escala do eixo y do gráfico
        self.fit_ax.set_title('Intensidade em função de Q') #define o título do gráfico
        self.fit_ax.set_ylabel('Intensidade') #define o texto do eixo y do gráfico
        self.fit_ax.set_xlabel('$Q [\u212B^{-1}]$') #u212B == Angstrom
        self.fit_ax.legend(loc='upper left', prop={'size':8})
        self.fit_sim_figure.show()

    def fit_save_data_button_clicked(self):
        '''Essa função salva os dados do fit. Salva o fit report, a curva fitada e o resíduo.'''
        #TODO Alterar essa função para que ela salve esses arquivos da pasta dos dados. Atualmente ela esta os salvando na pasta do programa.
        #index = self.fit_tabs.currentIndex()
        initial = 0
        end = 0

        ###As variáveis "inital" e "end" são usada pois os valores do resisduo são retornados para todos os arquivos juntos num único array.
        ###Assim, essas variáveis são usadas para controlar o ínicio e fim de valores de resíduo referente à cada arquivo. 

        for i in range(len(self.ff_file)): #percorre todos o arquivos de dados
            end += len(self.ff_file[i].q)
            if self.file_fit_cb[i].isChecked(): #verifica se o dado foi selecionado para o fit
                data = np.column_stack([self._line3[i].get_xdata(), self._line3[i].get_ydata()]) #salva os dados da curva ajustada na variável data
                path = str(self.ff_file[i].path.parent)+'\\'
                try:
                    np.savetxt(path+self.label_s[i]+'-fit-data.fit', data)
                    print('Fit data saved at '+path)
                    with open(path+self.label_s[i]+'-fit-report.out', 'w') as file:
                        file.write(lmf.fit_report(self.out))
                    print('Fit report saved at '+path)
                    np.savetxt(path+self.label_s[i]+'-residuo.res', np.column_stack([self.ff_file[i].q, self.out.residual[initial:end]]))
                    #salva os arquivos e printa no terminal o caminho para onde foram salvos
                except Exception as e:
                    print(e)
                    print("Couldn't save data")
            initial = end
                
    ################     FILE FUNCTIONS     ##################
    def addFile(self): 
        '''Esta função adicona os widgets necessários para importar arquivo'''
        self.file_final_layout = QtWidgets.QVBoxLayout()
        self.file_browse = QtWidgets.QPushButton('Search') #Botão de buscar arquivo       
        self.file_browse.setMaximumWidth(80) 
        self.file_browse.clicked.connect(self.open_file_dialog)
        self.file_final_layout.addWidget(self.file_browse)
        self.prepare_file_widgets()
        #self.add_file_tabs()
        self.file_final_layout.setAlignment(QtCore.Qt.AlignTop)
        self.file_remove_tab_button_group.buttonClicked.connect(self.remove_file_tabs)
        self.file_adj_button_group.buttonClicked.connect(self.file_ajd_button_clicked)        

    def prepare_file_widgets(self):
        '''Essa função cria as variavéis necessárias para a cração dos widgets de opções dos arquivos importados.'''
        self.file_adj = None
        self.file_layout = []
        self.file_groupbox = []
        self.filename_edit = []
        self.file_gb_layout = []
        self.file_remove_tab_button = []
        self.ff_file = []
        self.file_index = -1
        self.file_index_array = []
        self._line2 = []
        self._line2_error = []
        self._line3 = []
        self.file_fit_cb = []
        self.file_fit_label = []
        self.file_show_data_cb = []
        self.file_show_data_label = []
        self.file_q_unit_label = []
        self.file_q_unit_combobox = []
        self.file_q_unit_options = {'angstrom':'\u212b^-1', 'nano':'\u03b7m^-1', 'micro':'\u03bcm^-1'}

        self.file_remove_tab_button_group = QtWidgets.QButtonGroup()

        self.file_adj_button = []
        self.file_adj_button_group = QtWidgets.QButtonGroup()

        self.label_s = []
        self.color_s = []
        self.marker_s = []       # linestyle or marker
        self.s_s = []       # linewidth or size
        self.multiplier_s = []     # multiplier so multiple curves can be shifted whitout changing the intensity
        self.indmin_s, self.indmax_s, self.indmin_f, self.indmax_f = [], [], [], []     # minimum and maximum indexes for plotting and fitting
        self.valuemin_s, self.valuemax_s, self.valuemin_f, self.valuemax_f = [], [], [], []
        self.NameScatter_label = {}
    
    def add_file_tabs(self):
        '''Essa função cria os widgets da interface referente aos arquivos de dados importados.'''
        self.file_remove_tab_button.append(QtWidgets.QPushButton()) #botão para remover o arquivo
        self.file_remove_tab_button[-1].setIcon(qta.icon('mdi6.note-remove')) #seta o icone no botão
        self.file_remove_tab_button[-1].setToolTip('Remove File') #adiciona tooltip
        self.file_gb_layout.append(QtWidgets.QVBoxLayout()) 

        self.file_index += 1
        self.file_index_array.append(self.file_index)
        
        self.file_layout.append(QtWidgets.QGridLayout())
        self.file_groupbox.append(QtWidgets.QGroupBox('Arquivo'))
        self.filename_edit.append(QtWidgets.QLineEdit()) #caixa de texto com o caminho e nome do arquivo

        self.file_fit_cb.append(QtWidgets.QCheckBox()) #checkbox para selecionar se o arquivo será fitado ou não
        self.file_fit_label.append(QtWidgets.QLabel("Fit data: "))
        self.file_fit_label[-1].setMinimumWidth(40)
        self.file_fit_label[-1].setMaximumWidth(40)
        self.file_fit_cb[-1].setToolTip('Fit On/Off')
        self.file_fit_cb[-1].setChecked(True)
        self.file_fit_cb[-1].setMinimumWidth(20)
        self.file_fit_cb[-1].setMaximumWidth(20)

        self.file_show_data_cb.append(QtWidgets.QCheckBox()) #checkbox para exibir o dado no gráfico ou não
        self.file_show_data_label.append(QtWidgets.QLabel("Show data: "))
        self.file_show_data_label[-1].setMinimumWidth(53)
        self.file_show_data_label[-1].setMaximumWidth(53)
        self.file_show_data_cb[-1].setToolTip('Show/Hide data')
        self.file_show_data_cb[-1].setChecked(True)
        self.file_show_data_cb[-1].setMinimumWidth(20)
        self.file_show_data_cb[-1].setMaximumWidth(20)

        self.file_q_unit_label.append(QtWidgets.QLabel("Q unit: "))
        self.file_q_unit_label[-1].setMinimumWidth(33)
        self.file_q_unit_label[-1].setMaximumWidth(33)
        self.file_q_unit_combobox.append(QtWidgets.QComboBox())
        self.file_q_unit_combobox[-1].addItems(list(self.file_q_unit_options.values()))
        self.file_q_unit_combobox[-1].setMinimumWidth(55)
        self.file_q_unit_combobox[-1].setMaximumWidth(55)

        self.file_adj_button.append(QtWidgets.QPushButton()) #botão que abre janela de configurações do arquivo
        self.file_adj_button[-1].setIcon(self.dist_adj_icon) #seta icone para o botão
        self.file_adj_button[-1].setToolTip('Data/Fit options') 
        self.file_adj_button[-1].setMaximumWidth(40)
        
        self.filename_edit[-1].setMinimumWidth(150)
        self.filename_edit[-1].setMaximumWidth(150)

        #adiciona os widgets ao layout
        self.file_layout[-1].addWidget(self.file_show_data_label[-1], 0, 0)
        self.file_layout[-1].addWidget(self.file_show_data_cb[-1], 0, 1)
        self.file_layout[-1].addWidget(self.file_fit_label[-1], 0, 2)
        self.file_layout[-1].addWidget(self.file_fit_cb[-1], 0, 3)
        self.file_layout[-1].addWidget(self.file_q_unit_label[-1], 0, 4)
        self.file_layout[-1].addWidget(self.file_q_unit_combobox[-1], 0, 5)
        self.file_layout[-1].addWidget(self.file_adj_button[-1], 0, 6)
        self.file_layout[-1].addWidget(self.filename_edit[-1], 0, 7)
        self.file_layout[-1].addWidget(self.file_remove_tab_button[-1], 0, 8)

        #seta as funções a serem chamadas quando os botões são clicados
        self.file_remove_tab_button_group.addButton(self.file_remove_tab_button[-1], len(self.file_gb_layout)-1)
        self.file_adj_button_group.addButton(self.file_adj_button[-1], len(self.file_gb_layout)-1)

        self.file_q_unit_combobox[-1].currentIndexChanged.connect(lambda: self.file_q_unit_combobox_changed(len(self.file_q_unit_combobox)-1))

        self.file_show_data_cb[-1].stateChanged.connect(self.file_show_data_cb_changed)
        #print(self.file_remove_tab_button_group.id(self.file_remove_tab_button[-1]))
        self.file_remove_tab_button[-1].setMaximumWidth(50)

        #self.file_layout[-1].setAlignment(QtCore.Qt.AlignLeft)
        #self.file_layout[-1].setVerticalSpacing(2)
        #self.file_layout[-1].addStretch()

        self.file_groupbox[-1].setLayout(self.file_layout[-1])
        self.file_groupbox[-1].setMaximumHeight(80)

        self.file_gb_layout[-1].addWidget(self.file_groupbox[-1])
        
        self.ff_file.append(False)

        self.marker_s.append('o')
        self.s_s.append(5)
        self.multiplier_s.append(1)   
        self.file_final_layout.addWidget(self.file_groupbox[-1])

        self.any_file_added = True #seta como 'True' variável de controle para saber se algum arquivo de dados foi importado

    def file_q_unit_combobox_changed(self, file_index):
        unit = {'\u212b^-1':10e-10, '\u03b7m^-1':10e-9, '\u03bcm^-1':10e-6}
        self.ff_file[file_index].q_unit = self.file_q_unit_combobox[file_index].currentText()
        self.ff_file[file_index].q = self.ff_file[file_index].q_original*(unit[self.ff_file[file_index].q_unit]/self.graph_q_unit)
        self._update_file_line(file_index)
    
    def file_ajd_button_clicked(self, button):
        '''Essa função abre a jenela de configuração dos dados importados'''
        if self.file_adj == None:
            self.file_adj = PopUpDataOpt(self, self.file_adj_button_group.id(button))
            self.file_adj.show()
        elif self.file_adj.isVisible():
            self.file_adj.setWindowState(QtCore.Qt.WindowNoState)
        else:
            self.file_adj = PopUpDataOpt(self, self.file_adj_button_group.id(button))
            self.file_adj.show()
    
    def add_relative_file_cb_options(self):
        '''Essa função adicona os arquivos importados na combobox de seleção de arquivo relativo nas abas de fatores de forma'''
        i = self.ff_tabs.currentIndex()
        self.relative_file[i].addItem('None')
        for j in range(len(self.file_gb_layout)):
            self.relative_file[i].addItem(self.label_s[j])
        self.relative_file[i].setCurrentIndex(0)
    
    def upd_added_relative_file_cb_options(self):
        '''Essa função atualiza os arquivos importados na combobox de seleção de arquivo relativo nas abas de fatores de forma'''
        for i in range(len(self.tabs)):
            self.relative_file[i].addItem(self.label_s[-1])
    
    def upd_removed_relative_file_cb_options(self, index=0):
        '''Essa função remove os arquivos importados na combobox de seleção de arquivo relativo nas abas de fatores de forma'''
        for i in self.relative_file:
            #currentIndex = self.relative_file[i].currentIndex()
            i.removeItem(index)
            #for j in range(index, len(self.file_gb_layout)+1):
                #i.setItemText(j, self.label_s[j-1])
    
    def upd_label_changed_file_cb_options(self, option):
        for i in self.relative_file:
            i.setItemText(option+1, self.label_s[option])

    def file_show_data_cb_changed(self):
        '''Essa função atualiza a visibilidade das curvas importadas baseado na checkbox de visibilidade de cada arquivo.'''
        #TODO do jeito que está agora essa função percorre por todos os checboxs para verificar se estão selecionados ou não.
        for i in range(len(self.file_show_data_cb)):
            if self.file_show_data_cb[i].isChecked():
                self._line2[i].set_visible(True)
                #self._line2_error[i].set_visible(True)
            else:
                self._line2[i].set_visible(False)
                #self._line2_error[i].set_visible(False)
        self._line2[i].figure.canvas.draw() 

    def remove_file_tabs(self, button):
        '''Essa função remove um determinado arquivo importado.
        A variável "button" informa qual botão de remoção foi clicado para remoção do arquivo referente aquele botão'''
        tab = self.file_remove_tab_button_group.id(button) #salva o id do botão para saber qual arquivo remover

        self.file_final_layout.removeWidget(self.file_groupbox[tab]) 
        self.file_groupbox[tab].deleteLater()

        self.remove_file_widgets(tab)

        #for i in range(len(self.file_layout)):
            #print(self.file_layout[i].index(), self.file_remove_tab_button_group.id(self.file_remove_tab_button[i]))

        for i in range(tab, len(self.file_groupbox)): #percorre os arquivos que foram importados depois do arquivo que foi removido
            #troca o id desses arquivos considerando a remoção
            self.file_remove_tab_button_group.setId(self.file_remove_tab_button[i], i)
            self.file_adj_button_group.setId(self.file_adj_button[i], i)

        self.upd_removed_relative_file_cb_options(tab+1) #tira esse arquivo das opções de arquivos a serem selecionados nas abas de fator de forma
        self.ff_file.pop(tab)
        #self.remove_fit_tab(tab)

        self.create_lmf_params()

        self.valuechange('scale')

    def remove_file_widgets(self, index):
        '''Essa função remove os widgets de controle de um arquivo de dados'''
        self.file_layout.pop(index)
        self.file_groupbox.pop(index)
        self.filename_edit.pop(index)
        self.file_gb_layout.pop(index)
        self.file_remove_tab_button.pop(index)
        self._line2.pop(index).remove()
        self._line2_error.pop(index).remove()
        self._line3.pop(index).remove()
        self._line.pop(index).remove()
        if len(self.file_layout) == 0:
            self.any_file_added = False
        self.valuemin_s.pop(index)
        self.valuemax_s.pop(index)
        self.valuemin_f.pop(index)
        self.valuemax_f.pop(index)
        self.indmin_s.pop(index)
        self.indmax_s.pop(index)
        self.indmin_f.pop(index)
        self.indmax_f.pop(index)
        self.marker_s.pop(index)
        self.s_s.pop(index)
        self.multiplier_s.pop(index)
        self.label_s.pop(index)
        self.color_s.pop(index)
        self.file_fit_cb.pop(index)
        self.file_fit_label.pop(index)
        self.file_show_data_cb.pop(index)
        self.file_show_data_label.pop(index)
        self.file_adj_button.pop(index)
        self.file_q_unit_combobox.pop(index)
        self.file_q_unit_label.pop(index)

    def relative_file_changed(self, index = 0):
        '''Essa função é executada quando o arquivo relativo de alguma das abas de modelo de fator de forma e trocado.
        A variável "index" informa o número da aba onde o arquivo relativo foi trocado.'''
        #print(self.relative_file[index].currentIndex())
        if self.any_file_added and self.relative_file[index].currentIndex()-1>=0:
            #print('len(ff_file)-1={}, index={}'.format(len(self.ff_file)-1,self.relative_file[index].currentIndex()-1))
            if self.ff_file[self.relative_file[index].currentIndex()-1]:
                self.ffs[index].setq(self.ff_file[self.relative_file[index].currentIndex()-1].q)
                self.ffs[index].ff()
                if self.sf_added:
                    self.sf[index].setq(self.ff_file[self.relative_file[index].currentIndex()-1].q)
                    #self.sf[index].sf()
                self.valuechange('scale')

    def _update_file_line(self, scatter): 
        '''Essa função atualiza a formatação da curva dos dados importados.
        A variável "scatter" informa qual curva deve ser atualizada '''
        #####################################################
        _line2 = self.ax.scatter(self.ff_file[scatter].q[self.indmin_s[scatter]:self.indmax_s[scatter]], 
                                 self.multiplier_s[scatter]*self.ff_file[scatter].intensidade[self.indmin_s[scatter]:self.indmax_s[scatter]],
                                 marker=self.marker_s[scatter], label=self.label_s[scatter], s=self.s_s[scatter], color=self.color_s[scatter])
        _line2_error = self.ax.errorbar(self.ff_file[scatter].q[self.indmin_s[scatter]:self.indmax_s[scatter]],
                                        self.multiplier_s[scatter]*self.ff_file[scatter].intensidade[self.indmin_s[scatter]:self.indmax_s[scatter]],
                                        self.ff_file[-1].yerr[self.indmin_s[scatter]:self.indmax_s[scatter]], ecolor=self.color_s[scatter], fmt='none')
        self._line2[scatter].remove()
        self._line2[scatter] = _line2
        self._line2_error[scatter].remove()
        self._line2_error[scatter] = _line2_error
        self._line2[scatter].figure.canvas.draw() 
    
    def open_file_dialog(self): 
        '''Esta função importa os arquivos'''
        filename, ok = QtWidgets.QFileDialog.getOpenFileName(self, "Selecione um arquivo") #Cria a janela para seleção do arquivo.
        #Quando o arquivo é selecionado seu caminho fica salvo na variável "filename". Se o botão "Cancelar" é apartado dentro da janela de seleção a variável
        #"filename será nula"
        if filename: #se algum arquivo foi selecionado:
            self.add_file_tabs() #Adiciona na interface o widgets para configuração dos arquivos importados
            path = Path(filename)
            self.label_s.append(path.stem)
            self.filename_edit[-1].setText(str(path))
            #x, y = np.transpose(np.loadtxt(filename))
            #x, y = np.transpose(np.loadtxt(filename))
            for i in range(10):
                try:
                    print(i)
                    data = np.transpose(np.loadtxt(filename, skiprows=i))
                    break
                except:
                    pass
            print(len(data), data.shape)
            x = data[0]
            y = data[1]
            if len(data) == 4:
                y_error = data[3]
                x_error = data[2]
            elif len(data) == 3:
                y_error = data[2]
            elif len(data) == 2:
                y_error = np.zeros(len(data[1]))
            self.indmin_s.append(0)
            self.indmax_s.append(len(y))
            self.indmin_f.append(0)
            self.indmax_f.append(len(y))
            self.valuemin_s.append(min(x))
            self.valuemax_s.append(max(x))
            self.valuemin_f.append(min(x))
            self.valuemax_f.append(max(x))
            #self.color_s.append('#'+"%06x" % random.randint(0, 0xFFFFFF))
            self.ff_file[-1] = f.FileData(x, y, y_error)
            self.ff_file[-1].path = path
            '''for i in range(len(self.tabs)):
                if self.relative_file[i].currentIndex()-1 == self.file_tabs.currentIndex():
                    self.ffs[i].setq(x)               
                    self.ffs[i].ff() #calcula o fator de forma com o novo vetor q
                if i != len(self.tabs)-1:
                    self.valuechange('scale', i, False)
                else:
                    self.valuechange('scale', i)'''
            _line2, = self.ax.plot(self.ff_file[-1].q, self.ff_file[-1].intensidade, '.', label=Path(filename).stem, markersize=5)
            self.color_s.append(_line2.get_color())
            _line2_error = self.ax.errorbar(self.ff_file[-1].q, self.ff_file[-1].intensidade, self.ff_file[-1].yerr, ecolor=self.color_s[-1], fmt='none')
            _line3, = self.ax.plot(self.ff_file[-1].q, self.ff_file[-1].intensidade, label=Path(filename).stem+'fit2', color=self.colors.darken(self.color_s[-1]))
            ###_line3, = self.ax.plot(self.ff_file[-1].q, self.ff_file[-1].intensidade, color=self.colors.darken(self.color_s[-1]))
            self._line2.append(_line2)
            self._line2_error.append(_line2_error)
            self._line3.append(_line3)
            self._line3[-1].set_visible(False)
            self.NameScatter_label.update({len(self._line2)-1:QtWidgets.QLabel()})
            #self.ax.set_xlim([10**(np.floor(np.log10(min(x)))-0.5), 10**(np.floor(np.log10(max(self.ffs[self.ff_tabs.currentIndex()].q)))+0.5)])
            #self.updYLim()
            #for i in self._line2: i.figure.canvas.draw()    
            if len(self.file_layout)>0: #and not(self.first_file_execution):
                _line, = self.ax.plot(1, 1, '-', label=self.label_s[-1]+' fit', color=self.color_s[-1], markersize=5) #salva a curva do gráfico na variável line
                ###_line, = self.ax.plot(1, 1, '-', color=self.color_s[-1], markersize=5)
                self._line.append(_line)
            #print(self.file_layout.index(self.file_layout[-1]), self.file_remove_tab_button_group.id(self.file_remove_tab_button[-1]))
            self.first_file_execution = False
            self.upd_added_relative_file_cb_options()
            self.ax.legend(loc='upper left', prop={'size':8})
            #for i in self._line2: i.figure.canvas.draw()
            #for i in self._line3: i.figure.canvas.draw()
            self._line2[-1].figure.canvas.draw()
            #self.raioDeGiro() #calcula o raio de giro 
            #self.add_fit_tab()

    ####################   DISTRIBUITION FUNCTIONS    ####################
    def prepare_dist_widgets(self):
        self.dist_groupbox = []
        self.dist_layout = []
        self.dist_text = []
        self.dist_param_text = []
        self.dist_adj_button = []
        self.dist_box = []
        self.param_box = []
        self.dist_cb = []

    def add_dist(self): 
        '''Esta função adiciona a opção de aplicar distribuição aos parâmetros'''
        self.prepare_dist_widgets()
        self.add_dist_tab()
        self.dist_widgets_added[0] = True               

    def add_dist_tab(self):
        index = self.ff_tabs.currentIndex()

        self.dist_groupbox.append(QtWidgets.QGroupBox('Distribution'))
        self.dist_layout.append(QtWidgets.QGridLayout())
        self.dist_text.append(QtWidgets.QLabel('Distribution: '))
        self.dist_param_text.append(QtWidgets.QLabel('Parameter: '))
        self.dist_adj_button.append(QtWidgets.QPushButton())
        self.dist_adj_icon = qta.icon('ei.adjust-alt')
        #print(len(self.dist_adj_button)-1, index)
        self.dist_adj_button[index].setIcon(self.dist_adj_icon) #botão para abrir a jenala de configurações avançadas da distribuição a ser calculada
        self.dist_box.append(QtWidgets.QComboBox())
        self.param_box.append(QtWidgets.QComboBox())
        self.dist_cb.append(QtWidgets.QCheckBox())
        self.dist_box[index].setMaximumWidth(100)
        self.param_box[index].setMaximumWidth(100)
        self.dist_adj_button[index].setMaximumWidth(40)
        self.dist_box[index].addItems(['Gaussian', 'Log-Normal']) #opções de distribuição
        self.param_box[index].addItems(self.ffs[index].dist_params[0]) #opções de parâmetros para aplicar a distribuição

        self.dist_added = True #variável para controlar se a distribuição está sendo aplicada ou não

        self.dist_box[index].currentIndexChanged.connect(self.dist_box_changed)
        self.param_box[index].currentIndexChanged.connect(self.param_box_changed)
        self.dist_cb[index].stateChanged.connect(self.dist_cb_changed)
        self.dist_adj_button[index].clicked.connect(self.dist_adj_button_clicked)

        self.dist_layout[index].addWidget(self.dist_text[index], 1, 0)
        self.dist_layout[index].addWidget(self.dist_box[index], 1, 1)
        self.dist_layout[index].addWidget(self.dist_param_text[index], 2, 0)
        self.dist_layout[index].addWidget(self.param_box[index], 2, 1)
        self.dist_layout[index].addWidget(self.dist_cb[index], 0, 0)
        self.dist_layout[index].addWidget(self.dist_adj_button[index], 0, 1)

        self.addparam('dist', 0, 1, 0, 'Sigma', self.dist_layout[index]) #adiciona o slider do sigma da distribuição
        #self.change_param_position('dist', 3, 1, 1, 1, 'top', self.dist_layout)
        self.dist_cb_changed()

        self.dist_groupbox[index].setLayout(self.dist_layout[index])

        self.dist_groupbox[index].setMaximumHeight(210)

        self.gb_layout[index].addWidget(self.dist_groupbox[index],1)
        
    def remove_dist_widgets(self):
        index = self.ff_tabs.currentIndex()

        self.dist_groupbox.pop(index)
        self.dist_layout.pop(index)
        self.dist_text.pop(index)
        self.dist_param_text.pop(index)
        self.dist_adj_button.pop(index)
        self.dist_box.pop(index)
        self.param_box.pop(index)
        self.dist_cb.pop(index)
    
    def param_box_changed(self):
        if self.param_box[self.ff_tabs.currentIndex()].currentIndex() >= 0:
            self.valuechange('scale')

    def dist_cb_changed(self):
        '''Esta função habilita ou desabilita a opção de adicionar a distribuição'''
        index = self.ff_tabs.currentIndex()
        var = self.dist_cb[index].isChecked()

        self.dist_box[index].setEnabled(var)
        self.param_box[index].setEnabled(var)
        self.slider['Model'][index]['dist'].setEnabled(var)
        self.text_edit_max['Model'][index]['dist'].setEnabled(var)
        self.text_edit_min['Model'][index]['dist'].setEnabled(var)
        self.text_edit_value['Model'][index]['dist'].setEnabled(var)
        self.dist_adj_button[index].setEnabled(var)
        self.param_vary_cb['Model'][index]['dist'].setEnabled(var)
        self.param_expr_textbox['Model'][index]['dist'].setEnabled(var)
        self.valuechange('dist')

    def dist_box_changed(self):
        self.valuechange('dist')

    def param_box_changed(self):
        self.valuechange('dist')

    def dist_adj_button_clicked(self):
        if self.adj == None:
            self.adj = DistAdj(self)
            self.adj.show()
        elif self.adj.isVisible():
            self.adj.setWindowState(QtCore.Qt.WindowNoState)
        else:
            self.adj = DistAdj(self)
            self.adj.show()

    ######################     STRUCTURE FACTORS FUNCTIONS      ###############################
    def prepare_sf_widgets(self):
        self.sf_layout = []
        self.sf_groupbox = []
        self.sf = []
        self.sf_text = []
        self.sf_param_text = []
        self.sf_aprox_text = []
        #self.sf_Ref_text = []
        self.sf_box = []
        self.sf_aprox_box = []
        self.sf_param_box = []
        self.sf_cb = []

    def add_sf(self): 
        '''Esta função adiciona a opção de aplicar o fator de estrutura'''
        self.prepare_sf_widgets()
        self.add_sf_tab()
        self.sf_widgets_added[0] = True

    def add_sf_tab(self):
        index = self.ff_tabs.currentIndex()

        self.sf_layout.append(QtWidgets.QGridLayout())
        self.sf_groupbox.append(QtWidgets.QGroupBox('Structure Factor'))
        self.sf.append(s.HardSphere(50))
        self.sf[index].setq(self.ffs[index].q)
        self.sf_text.append(QtWidgets.QLabel('Structure Factor: '))
        self.sf_param_text.append(QtWidgets.QLabel('Parameter: '))
        self.sf_aprox_text.append(QtWidgets.QLabel('Approximation: '))
        #self.sf_Ref_text.append(QtWidgets.QLabel('Effective Radius: '))
        #self.sf_Ref_textbox[index] = QtWidgets.QLineEdit(list(self.ffs[index].func_pars.keys())[0])
        self.param_expr_textbox['Model'][index]['r_ef'] = list(self.ffs[index].func_pars.keys())[0]+str(index)
        #self.sf_adj_button = QtWidgets.QPushButton()
        #self.sf_adj_icon = qta.icon('ei.adjust-alt')
        #self.sf_adj_button.setIcon(self.dist_adj_icon) #botão para abrir a jenala de configurações avançadas da distribuição a ser calculada
        self.sf_box.append(QtWidgets.QComboBox())
        self.sf_aprox_box.append(QtWidgets.QComboBox())
        self.sf_param_box.append(QtWidgets.QComboBox())
        self.sf_cb.append(QtWidgets.QCheckBox())
        self.sf_box[index].setMaximumWidth(100)
        self.sf_param_box[index].setMaximumWidth(100)
        self.sf_aprox_box[index].setMaximumWidth(100)
        #self.sf_Ref_textbox[index].setMaximumWidth(100)
        #self.sf_adj_button.setMaximumWidth(40)
        self.sf_box[index].addItems(s.all) #opções de distribuição
        self.sf_param_box[index].addItems(self.sf[index].methods) #opções de parâmetros para aplicar a distribuição
        self.sf_aprox_box[index].addItems(s.aprox) #opções de parâmetros para aplicar a distribuição

        self.sf_added = True #variável para controlar se a distribuição está sendo aplicada ou não

        self.sf_box[index].currentIndexChanged.connect(self.sf_box_changed)
        self.sf_param_box[index].currentIndexChanged.connect(self.sf_param_box_changed)
        self.sf_cb[index].stateChanged.connect(self.sf_cb_changed)
        #self.sf_Ref_textbox[index].returnPressed.connect(self.sf_Ref_textbox_pressed)
        self.sf_aprox_box[index].currentIndexChanged.connect(self.sf_aprox_changed)
        #self.sf_adj_button.clicked.connect(self.sf_adj_button_clicked)

        self.sf_layout[index].addWidget(self.sf_aprox_text[index], 1, 0)
        self.sf_layout[index].addWidget(self.sf_aprox_box[index], 1, 1)
        self.sf_layout[index].addWidget(self.sf_text[index], 2, 0)
        self.sf_layout[index].addWidget(self.sf_box[index], 2, 1)
        self.sf_layout[index].addWidget(self.sf_param_text[index], 3, 0)
        self.sf_layout[index].addWidget(self.sf_param_box[index], 3, 1)
        self.sf_layout[index].addWidget(self.sf_cb[index], 0, 0)
        #self.sf_layout[index].addWidget(self.sf_Ref_text[index], 4, 0)
        #self.sf_layout[index].addWidget(self.sf_Ref_textbox[index], 4, 1)

        #self.layout.addWidget(self.sf_adj_button, self.position2+1, 1)
        position = 5
        for i in self.sf[index].func_pars:
            self.addparam(i, self.sf[index].func_pars[i]['min'], self.sf[index].func_pars[i]['max'],
                            self.sf[index].func_pars[i]['value'], self.sf[index].func_pars[i]['symbol'], 
                            self.sf_layout[index])
            #elf.change_param_position(i, position, 1, 1, 1, 'top')
            position+=1
        self.sf_cb_changed()

        self.sf_groupbox[index].setLayout(self.sf_layout[index])
        self.sf_groupbox[index].setMaximumHeight(240)
        self.gb_layout[index].addWidget(self.sf_groupbox[index], 4)
        
    def remove_sf_widgets(self):
        index = self.ff_tabs.currentIndex()

        self.sf_layout.pop(index)
        self.sf_groupbox.pop(index)
        self.sf.pop(index)
        self.sf_text.pop(index)
        self.sf_param_text.pop(index)
        self.sf_aprox_text.pop(index)
        #self.sf_Ref_text.pop(index)
        self.sf_box.pop(index)
        self.sf_aprox_box.pop(index)
        self.sf_param_box.pop(index)
        self.sf_cb.pop(index)
        #self.sf_Ref_textbox.pop(index)
    
    def sf_cb_changed(self):
        '''Esta função habilita ou desabilita a opção de adicionar o fator de estrutura'''
        index = self.ff_tabs.currentIndex()
        var = self.sf_cb[index].isChecked()

        self.sf_box[index].setEnabled(var)
        self.sf_param_box[index].setEnabled(var)
        self.sf_aprox_box[index].setEnabled(var)
        #self.sf_Ref_textbox[index].setEnabled(var)
        for i in self.sf[index].func_pars:
            self.slider['Model'][index][i].setEnabled(var)
            self.text_edit_max['Model'][index][i].setEnabled(var)
            self.text_edit_min['Model'][index][i].setEnabled(var)
            self.text_edit_value['Model'][index][i].setEnabled(var)
            self.param_vary_cb['Model'][index][i].setEnabled(var)
            self.param_expr_textbox['Model'][index][i].setEnabled(var)
        #self.sf_adj_button.setEnabled(True)
        self.valuechange('fp')

    def sf_box_changed(self):
        ''''''
        self.valuechange('fp')

    def sf_param_box_changed(self):
        self.valuechange('fp')

    def sf_Ref_textbox_pressed(self):
        self.valuechange('fp')
    
    def sf_aprox_changed(self):
        self.valuechange('fp')

    def sf_ref_expr(self, index = None, funcpars=False, funcpars_sf = False):
        #! Essa função faz o papel do expression do lmfit
        #! Foi criada pois antes eu não usava o lmfit para salvar os valores dos parametros
        #! Agora é necessário reavaliar se ela ainda faz algum sentido
        if index == None:
            index = self.ff_tabs.currentIndex()
        
        if funcpars:
            for i in funcpars:
                #print(i)
                exec(i+"= funcpars[i]['value']")
            if funcpars_sf:
                for j in funcpars_sf:
                    exec(j+"= funcpars_sf[j]['value']")
        else:
            for i in self.slider['Model'][index]:
                #print(i)
                exec(i+'= self.slider["Model"][index][i].value()')
        #print(len(self.expr), len(self.sf_Ref_textbox), index)
        self.expr[index] = split(r'[^\w]', self.sf_Ref_textbox[index].text())
        if self.dist_added:
            if self.dist_cb[index].isChecked():
                if self.slider['Model'][index]['dist'].value() != 0 and (self.ffs[index].p in self.expr[index]) and self.sf_aprox_box[index].currentIndex()==1:
                    if self.dist_cb[index].isChecked():
                        exec(self.ffs[index].p+'= self.ffs[index].dist.x')
                        r_ef = eval(self.sf_Ref_textbox[index].text())
                    else:
                        r_ef = eval(self.sf_Ref_textbox[index].text())
                else:
                    r_ef = eval(self.sf_Ref_textbox[index].text())
            else:
                r_ef = eval(self.sf_Ref_textbox[index].text())
        else:
            r_ef = eval(self.sf_Ref_textbox[index].text())
        return r_ef

    ####################        PEAKS        #########################
    def peak_prepare_widgets(self):
        self.peaks_tabs = QtWidgets.QTabWidget()
        self.peaks_tabs_array = []
        self.peaks_tabs_layouts = []
        self.peaks_final_layout = QtWidgets.QVBoxLayout()
        self.peaks_structure_combobox = []
        self.peaks_q_slider = []
        self.peaks_activate_cb = []
        self.peaks_structure = []
        self._line4 = []
        self.peaks_line_color = []
        self.peaks_previous_point_button = []
        self.peaks_next_point_button = []
        self.peaks_final_layout.addWidget(self.peaks_tabs)
        self.adj_tabs_array[3].setLayout(self.peaks_final_layout)

    def peak_add_tab(self):
        self.peaks_tabs_layouts.append(QtWidgets.QGridLayout())
        self.peaks_tabs_array.append(QtWidgets.QWidget())
        self.peaks_tabs_array[-1].setLayout(self.peaks_tabs_layouts[-1])
        self.peaks_tabs_layouts[-1].setAlignment(QtCore.Qt.AlignTop)
        self.peaks_tabs.addTab(self.peaks_tabs_array[-1], 'Tab {}'.format(len(self.peaks_tabs_layouts)))
        self.peaks_structure_combobox.append(QtWidgets.QComboBox())
        self.peaks_structure_combobox[-1].addItems(p.structures)
        self.peaks_activate_cb.append(QtWidgets.QCheckBox())
        self.peaks_previous_point_button.append(QtWidgets.QPushButton())
        self.peaks_previous_point_button[-1].setIcon(qta.icon('fa.arrow-left'))
        self.peaks_next_point_button.append(QtWidgets.QPushButton())
        self.peaks_next_point_button[-1].setIcon(qta.icon('fa.arrow-right'))
        self.peaks_activate_cb[-1].stateChanged.connect(self.peak_activate_cb_changed)
        self.peaks_structure_combobox[-1].currentIndexChanged.connect(self.peaks_sctructure_changed)
        self.peaks_tabs_layouts[-1].addWidget(self.peaks_activate_cb[-1], 0,0,1,1)
        self.peaks_tabs_layouts[-1].addWidget(self.peaks_structure_combobox[-1], 0, 1, 1, 2)
        self.peaks_tabs_layouts[-1].addWidget(self.peaks_previous_point_button[-1], 1, 1, 1, 1)
        self.peaks_tabs_layouts[-1].addWidget(self.peaks_next_point_button[-1], 1, 2, 1, 1)
        self.peaks_structure.append(self.which_structure())
        self.addNeededWidgets('Peaks')
        self.addparam('peak_q', 1e-4, 1, 1e-2, 'q', self.peaks_tabs_layouts[-1], 5, 'Peaks')
        self.change_param_position('peak_q', 1, 4, 1, 1, 'center', self.peaks_tabs_layouts[-1], 'Peaks')
        _line4 = self.ax.vlines([0.03], ymin=0, ymax=1e25)
        self._line4.append(_line4)
        self._line4[-1].set_visible(False)
        self.peaks_line_color.append(self._line4[-1].get_color())
        self.peak_activate_cb_changed()
        self.create_lmf_params()
        self.valuechange('peak_q', None, True, 'Peaks')

    '''def peaks_prev_or_next_button_clicked(self, button):
        if len(self.ff_file) == 0:
            if button == 'previous':
                self.slider[]
        difference_array = np.absolute(arr-x)
        index = difference_array.argmin()   ''' 
    
    def peak_activate_cb_changed(self):
        index = self.peaks_tabs.currentIndex()
        activate = self.peaks_activate_cb[index].isChecked()

        self.slider['Peaks'][index]['peak_q'].setEnabled(activate)
        self.text_edit_max['Peaks'][index]['peak_q'].setEnabled(activate)
        self.text_edit_min['Peaks'][index]['peak_q'].setEnabled(activate)
        self.text_edit_value['Peaks'][index]['peak_q'].setEnabled(activate)
        self.param_vary_cb['Peaks'][index]['peak_q'].setEnabled(activate)
        self.param_expr_textbox['Peaks'][index]['peak_q'].setEnabled(activate)
        self.peaks_structure_combobox[index].setEnabled(activate)
        self._line4[index].set_visible(activate)
        self._line4[index].figure.canvas.draw()

    def peaks_sctructure_changed(self):
        #todo fazer a função
        index = self.peaks_tabs.currentIndex()
        self.peaks_structure[index]=(self.which_structure(self.peaks_structure_combobox[index].currentIndex()))
        self.valuechange('peak_q', None, True, 'Peaks')
        #print(self.peaks_structure[index].name)

    def which_structure(self, structure=0):
        s = p.dic[structure]
        func = s()
        return func
    
    def upd_peaks(self, index):
        self.peaks_structure[index].lattice_parameter(self.lmf_params['peak_q'+'_comma_'+'Peaks'+'_comma_'+str(index)].value, self.peaks_structure[index].peaks_spacings)
        #print(self.lmf_params['peak_q'+str(index)].value)
        peaks_positions = self.peaks_structure[index].peak_position(self.peaks_structure[index].peaks_spacings)
        #print(peaks_positions)
        self._line4[index].set_segments([np.array([[x, self._line4[index].get_segments()[0][0,1]], [x, self._line4[index].get_segments()[0][1,1]]]) for x in peaks_positions])
        self._line4[index].set_color(self.peaks_line_color[index])
        self._line4[index].figure.canvas.draw()
        #_line4, = self.ax.vlines(self.peaks_structure[index].peak_position())#, ymin=min(self.ff_file[index].intensidade), ymax=max(self.ff_file[index].intensidade))
        #plot(self.ff_file[-1].q, self.ff_file[-1].intensidade, '.', label=Path(filename).stem, markersize=5)
                
    ####################     OTHER FUNCTIONS      ######################
    def flatten(self, l):
        return np.array([item for sublist in l for item in sublist])

    def format_e(self, n): 
        '''Esta função formata um valor para notação científica'''
        a = '%E' % n
        return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]
    
    def toolTip(self, name, tab): 
        '''Esta função adiciona uma legenda quando o mouse esta sobre a caixa de texto dos slider'''
        index = self.ff_tabs.currentIndex()
        
        self.text_edit_min[tab][index][name].setToolTip('Maximum Value')
        self.text_edit_max[tab][index][name].setToolTip('Minimum Value')
        self.text_edit_value[tab][index][name].setToolTip('Value')
        self.param_vary_cb[tab][index][name].setToolTip('Vary parameter on fit?')
        self.param_expr_textbox[tab][index][name].setToolTip('Expression')

    def alignment(self, name, alignment, tab): 
        '''Esta função aplica um alinhamento ao texto das caixas de texto do slider
        A variável 'name' é referente à chave do slider no dicinário de sliders'''   
        index = self.ff_tabs.currentIndex()

        if alignment == 'center':
            self.text_edit_min[tab][index][name].setAlignment(QtCore.Qt.AlignCenter)
            self.text_edit_max[tab][index][name].setAlignment(QtCore.Qt.AlignCenter)
            self.text_edit_value[tab][index][name].setAlignment(QtCore.Qt.AlignCenter)
            self.param_expr_textbox[tab][index][name].setAlignment(QtCore.Qt.AlignCenter)
        if alignment == 'right':
            self.text_edit_min[tab][index][name].setAlignment(QtCore.Qt.AlignRight)
            self.text_edit_max[tab][index][name].setAlignment(QtCore.Qt.AlignRight)
            self.text_edit_value[tab][index][name].setAlignment(QtCore.Qt.AlignRight)
            self.param_expr_textbox[tab][index][name].setAlignment(QtCore.Qt.AlignRight)
        if alignment == 'left':
            self.text_edit_min[tab][index][name].setAlignment(QtCore.Qt.AlignLeft)
            self.text_edit_max[tab][index][name].setAlignment(QtCore.Qt.AlignLeft)
            self.text_edit_value[tab][index][name].setAlignment(QtCore.Qt.AlignLeft)
            self.param_expr_textbox[tab][index][name].setAlignment(QtCore.Qt.AlignLeft)

    def change_param_position(self, name, row, column, lenght_row=1, lenght_column=1, alignment='center', layout = QtWidgets.QGridLayout(), tab = 'Model'): 
        '''Esta função muda a posição de um slider e suas caixa de texto dentro do layout da janela'''
        index = self.ff_tabs.currentIndex()

        if alignment == 'top':
            alignment = QtCore.Qt.AlignmentFlag.AlignTop
        elif alignment == 'center':
            alignment = QtCore.Qt.AlignmentFlag.AlignCenter
        elif alignment == 'bottom':
            alignment = QtCore.Qt.AlignmentFlag.AlignBottom
        elif alignment == 'right':
            alignment = QtCore.Qt.AlignmentFlag.AlignRight
        elif alignment == 'left':
            alignment = QtCore.Qt.AlignmentFlag.AlignLeft

        layout.addWidget(self.param_vary_cb[tab][index][name], row, column-1, alignment)
        layout.addWidget(self.param_text[tab][index][name], row, column, alignment)
        layout.addWidget(self.text_edit_min[tab][index][name], row, column+2, alignment)
        layout.addWidget(self.text_edit_max[tab][index][name], row, column+3, alignment)
        layout.addWidget(self.text_edit_value[tab][index][name], row, column+4, alignment) 
        layout.addWidget(self.slider[tab][index][name], row, column+1, lenght_row, lenght_column, alignment)
        layout.addWidget(self.param_expr_textbox[tab][index][name], row, column+5, alignment)

    def SS(self): #formatação visual
        return 'QLineEdit {' + 'border: 2px solid gray; border-radius: 10px; padding: 0 8px; selection-background-color: darkgray;}'
    
    def BS(self): #formatação visual
        return 'QPushButton {'+'border-style: outset;border-width: 2px;border-radius: 10px;border-color: gray;}'

class DistAdj(QtWidgets.QWidget):
    '''Classe da aba de ajuste da distribuição'''
    def __init__(self, main_class):
        super().__init__()
        self.apply_button = QtWidgets.QPushButton('Apply')
        self.cancel_button = QtWidgets.QPushButton('Cancel')
        self.plot_button = QtWidgets.QPushButton('Distribution')
        self.layout = QtWidgets.QGridLayout(self)
        self.npoints_text = QtWidgets.QLabel('N. of points: ')
        self.npoints_textbox = QtWidgets.QLineEdit()            
        self.min_x_text = QtWidgets.QLabel('Min x: ')
        self.max_x_text = QtWidgets.QLabel('Max x: ')
        self.min_sigma_text = QtWidgets.QLabel('\u03C3')
        self.max_sigma_text = QtWidgets.QLabel('\u03C3')
        self.apply_button.clicked.connect(lambda: self.apply_button_clicked(main_class))
        self.cancel_button.clicked.connect(lambda: self.close())
        self.plot_button.clicked.connect(lambda: self.plot_button_clicked(main_class))
        self.min_x_textbox = QtWidgets.QLineEdit()
        self.max_x_textbox = QtWidgets.QLineEdit()
        if main_class.ffs[main_class.ff_tabs.currentIndex()].dist:
            self.min_x_textbox.setText(str(main_class.ffs[main_class.ff_tabs.currentIndex()].dist.min_amp))
            self.npoints_textbox.setText(str(main_class.ffs[main_class.ff_tabs.currentIndex()].dist.n))
            self.max_x_textbox.setText(str(main_class.ffs[main_class.ff_tabs.currentIndex()].dist.max_amp))
        

        self.layout.addWidget(self.npoints_text, 0, 0)
        self.layout.addWidget(self.npoints_textbox, 0, 1)
        self.layout.addWidget(self.min_x_text, 1,0)
        self.layout.addWidget(self.min_x_textbox, 1,1)
        self.layout.addWidget(self.min_sigma_text, 1,2)
        #self.layout.addWidget(self.min_sigma_cb, 1,3)
        self.layout.addWidget(self.max_x_text, 2,0)
        self.layout.addWidget(self.max_x_textbox, 2,1)
        self.layout.addWidget(self.max_sigma_text, 2,2)
        #self.layout.addWidget(self.max_sigma_cb, 2,3)
        self.layout.addWidget(self.apply_button, 3, 4)
        self.layout.addWidget(self.cancel_button, 3, 3)
        self.layout.addWidget(self.plot_button, 3,2)

    def apply_button_clicked(self, main_class):
        if round(float(self.npoints_textbox.text())) > 100:
            self.warning_message(main_class)
        else:
            self.apply_changes(main_class)

    def apply_changes(self, main_class):
        main_class.min_amp = float(self.min_x_textbox.text())
        main_class.max_amp = float(self.max_x_textbox.text())
        main_class.n = round(float(self.npoints_textbox.text()))
        main_class.dist_adj_applied = True
        main_class.valuechange('dist')
    
    def plot_button_clicked(self, main_class):
        if main_class.slider['Model'][main_class.ff_tabs.currentIndex()]['dist'].value() == 0:
            self.error_message()
        elif main_class.dist_added:
            main_class.ffs[main_class.ff_tabs.currentIndex()].dist.plot()

    def error_message(self):
        self.msg = QtWidgets.QMessageBox()
        self.msg.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.msg.setIcon(QtWidgets.QMessageBox.Critical)
        self.msg.setText('Erro!')
        self.msg.setDetailedText("Para o cálculo da distribuição é necessário que sigma seja maior que zero.")
        self.msg.setWindowTitle('Erro')
        self.msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.msg.exec_()

    def warning_message(self, main_class):
        self.msg = QtWidgets.QMessageBox()
        self.msg.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.msg.setIcon(QtWidgets.QMessageBox.Warning)
        self.msg.setText('Para grande número de pontos o cálculo da distribuição e do fator de forma pode ser demorado.')
        self.msg.setWindowTitle('Aviso')
        self.msg.setStandardButtons(QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.Ok)
        self.msg.buttonClicked.connect(lambda: self.w_button_clicked(main_class))
        self.msg.exec_()

    def w_button_clicked(self, main_class):
        button = self.msg.clickedButton()
        sb = self.msg.standardButton(button)
        if sb == QtWidgets.QMessageBox.Ok:
            self.apply_changes(main_class)
        elif sb == QtWidgets.QMessageBox.Cancel:
            self.w_cancel_button(main_class)
    
    def w_cancel_button(self, main_class):
        #print(main_class.ff.dist)
        if main_class.ffs[0].dist:
            self.npoints_textbox.setText(str(main_class.ffs[0].dist.n))
        else:
            self.npoints_textbox.setText('')

class DoubleSlider(QtWidgets.QSlider): #slider
    # create our signal that we can connect to if necessary
    doubleValueChanged = QtCore.Signal(float) 
    
    def __init__(self, decimals=3, *args, **kargs):
        super(DoubleSlider, self).__init__( *args, **kargs)
        self._multi = 10 ** decimals
        self.valueChanged.connect(self.emitDoubleValueChanged)
        self.a = float(super(DoubleSlider, self).value())

    def emitDoubleValueChanged(self):
        value = float(super(DoubleSlider, self).value())/self._multi
        self.doubleValueChanged.emit(value) 

    def value(self):
        return float(super(DoubleSlider, self).value()) / self._multi 
    
    def minimum(self):
        return float(super(DoubleSlider, self).minimum()) / self._multi 
    
    def maximum(self):
        return float(super(DoubleSlider, self).maximum()) / self._multi 
    
    def func_value(self, a, i): #retorna o valor do slider aplicado a uma função que pode ser trocada
        return a

    def setMinimum(self, value):
        return super(DoubleSlider, self).setMinimum(int(value * self._multi)) 

    def setMaximum(self, value):
        return super(DoubleSlider, self).setMaximum(int(value * self._multi)) 

    def setSingleStep(self, value):
        return super(DoubleSlider, self).setSingleStep(int(value * self._multi)) 

    def singleStep(self):
        return float(super(DoubleSlider, self).singleStep()) / self._multi 

    def setValue(self, value):
        super(DoubleSlider, self).setValue(int(value * self._multi))

class PopUpDataOpt(QtWidgets.QWidget):
    def __init__(self, main_widget, scatter):
        QtWidgets.QWidget.__init__(self)
        
        self.layout = QtWidgets.QVBoxLayout(self)
        self.setWindowTitle("options")
        self.GroupBox = QtWidgets.QGroupBox('options')
        self.layout.addWidget(self.GroupBox)
        
        self.Groupbox_layout = QtWidgets.QGridLayout(self.GroupBox)
        self.Groupbox_layout.setAlignment(QtCore.Qt.AlignLeft)
        
        self.options = ['marker', 'size', 'label', 'multiplier', 'x min plot', 'x max plot', 'x min fit', 'x max fit']
        self.related_funcs = {'marker':main_widget.marker_s, 'size':main_widget.s_s, 'label':main_widget.label_s, 'multiplier':main_widget.multiplier_s,
                              'x min plot':main_widget.valuemin_s , 'x max plot':main_widget.valuemax_s, 'x min fit':main_widget.valuemin_f , 'x max fit':main_widget.valuemax_f}
        self.related_funcs2 = {'x min plot':main_widget.indmin_s , 'x max plot':main_widget.indmax_s, 'x min fit':main_widget.indmin_f , 'x max fit':main_widget.indmax_f}
        
        self.Labels = {}
        self.Entries = {}
        #print('relate_funcs: {}'.format(self.related_funcs))
        for j, i in enumerate(self.options):
            #print(j, i)
            self.Labels.update({i:QtWidgets.QLabel()})
            self.Labels[i].setText('{}'.format(i))
            self.Labels[i].setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
            
            self.Entries.update({i:QtWidgets.QLineEdit('{}'.format(self.related_funcs[i][scatter]))})
            self.Entries[i].setMaxLength(12)
            self.Entries[i].editingFinished.connect(lambda scat = scatter, opts = i: self.change(scat, opts))
            self.Entries[i].setFixedWidth(60)
            self.Entries[i].setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            self.Entries[i].setToolTip('enter new parameter value')

            self.Groupbox_layout.addWidget(self.Labels[i],j,0,1,1)
            self.Groupbox_layout.addWidget(self.Entries[i],j,1,1,1)
        
        
        self.main_widget = main_widget
    def change(self, scatter, name):
        self.markers = [".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "|", "_", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        if name == 'marker':
            if self.Entries[name].text() in self.markers:
                self.related_funcs[name][scatter] = self.Entries[name].text()
            else:
                print ('marker {} not implemented, sorry'.format(self.Entries[name].text()))
                self.Entries[name].setText(self.related_funcs[name][scatter])
        if name == 'size':
            try: 
                float_ = float(self.Entries[name].text())
            except:
                print ('size {} not a number, sorry'.format(self.Entries[name].text()))
                self.Entries[name].setText(self.related_funcs[name][scatter])
            if float_ > 40 or float_ < 1:
                print ('size {} not allowed, sorry'.format(self.Entries[name].text()))
                self.Entries[name].setText(str(self.related_funcs[name][scatter]))
            else:
                self.related_funcs[name][scatter] = float_
        if name == 'label':
            self.related_funcs[name][scatter] = self.Entries[name].text()
            self.main_widget.NameScatter_label[scatter].setText(self.Entries[name].text())
            self.main_widget._line2[scatter].set_label(self.Entries[name].text())
            self.main_widget._line[scatter].set_label(self.Entries[name].text()+' fit')
            self.main_widget._line3[scatter].set_label(self.Entries[name].text()+'fit2')
            self.main_widget.upd_label_changed_file_cb_options(scatter)
            self.main_widget.ax.legend(loc='upper left', prop={'size':8})
            self.main_widget._line[-1].figure.canvas.draw()
        if name == 'multiplier':
            try: 
                float_ = float(self.Entries[name].text())
                self.related_funcs[name][scatter] = float_
            except:
                print ('multiplier {} is not a number, sorry'.format(self.Entries[name].text()))
                self.Entries[name].setText(str(self.related_funcs[name][scatter]))
        if name in ['x min plot', 'x max plot', 'x min fit', 'x max fit']:
            try:
                float_ = float(self.Entries[name].text())
                self.related_funcs[name][scatter] = float_
                #print(min(self.main_widget.ff_file[scatter].q), float_)
                #print(np.searchsorted(min(self.main_widget.ff_file[scatter].intensidade), float_))
                self.related_funcs2[name][scatter] = np.searchsorted(self.main_widget.ff_file[scatter].q, float_)
            except Exception as e:
                print(e)
                print ('value {} is not a valid number, sorry'.format(self.Entries[name].text()))
                self.Entries[name].setText(str(self.related_funcs[name][scatter]))
            
        self.main_widget._update_file_line(scatter)
        self.main_widget.upd_graph()

class Colors():
    def __init__(self):
        pass
    def darken(self, color): #in format '#rrggbb'
        r = float(int(color[1:3],16))/255.
        g = float(int(color[3:5],16))/255.
        b = float(int(color[5:7],16))/255.
        h, s, v = colorsys.rgb_to_hsv(r,g,b)
        news = s 
        newv = v/2
        _newr, _newg, _newb = colorsys.hsv_to_rgb(h,news,newv)
        newr = hex(int(_newr*255))[2:4]
        newg = hex(int(_newg*255))[2:4]
        newb = hex(int(_newb*255))[2:4]
        return '#{:0>2}{:0>2}{:0>2}'.format(newr,newg,newb)
    def gray (self, color, value = 0.5):   #in format '#rrggbb' -> turning colors more gray. "value" should be between 0 (little gray) and 1 (gray)
        r = float(int(color[1:3],16))/255.
        g = float(int(color[3:5],16))/255.
        b = float(int(color[5:7],16))/255.
        h, s, v = colorsys.rgb_to_hsv(r,g,b)
        news = s - value
        if news < 0: news = 0
        newv = v
        _newr, _newg, _newb = colorsys.hsv_to_rgb(h,news,newv)
        newr = hex(int(_newr*255))[2:4]
        newg = hex(int(_newg*255))[2:4]
        newb = hex(int(_newb*255))[2:4]
        return '#{:0>2}{:0>2}{:0>2}'.format(newr,newg,newb)
    def lighten(self, color): #in format '#rrggbb'
        r = float(int(color[1:3],16))/255.
        g = float(int(color[3:5],16))/255.
        b = float(int(color[5:7],16))/255.
        h, s, v = colorsys.rgb_to_hsv(r,g,b)
        news = s/2. 
        newv = 0.5+v/2.
        _newr, _newg, _newb = colorsys.hsv_to_rgb(h,news,newv)
        newr = hex(int(_newr*255))[2:4]
        newg = hex(int(_newg*255))[2:4]
        newb = hex(int(_newb*255))[2:4]
        return '#{:0>2}{:0>2}{:0>2}'.format(newr,newg,newb)
    
class NewLmfParameters(lmf.Parameters):
    def __init__(self, main_class=None):
        super().__init__()
        self.main_class=main_class

    def update_constraints(self):
        #self.pretty_print()
        """Update all constrained parameters.

        This method ensures that dependencies are evaluated as needed.

        """
        requires_update = {name for name, par in self.items() if par._expr is not None}
        updated_tracker = set(requires_update)
        #print('1: ', requires_update)
        #print('2: ', updated_tracker)
        updated_tracker_ = set()
        for i in updated_tracker:
            i_ = i.replace('_comma_','').replace('Model', '').replace('Peaks','')
            updated_tracker_.add(i_)

        def _update_param(name):
            """Update a parameter value, including setting bounds.

            For a constrained parameter (one with an `expr` defined), this
            first updates (recursively) all parameters on which the
            parameter depends (using the 'deps' field).

            """
            #print('3: ', name)
            sucess = False
            for i in self.keys():
                #print(i.replace('_comma_','').replace('Model', '').replace('Peaks',''))
                if i.replace('_comma_','').replace('Model', '').replace('Peaks','') == name:
                    par = self[i]
                    sucess = True
            if sucess == False:    
                par = self.__getitem__(name)
            if par._expr_eval is None:
                par._expr_eval = self._asteval
            for dep in par._expr_deps:
                #print('4: ', dep)
                if dep in updated_tracker:
                    _update_param(dep)
                if dep in updated_tracker_:
                    _update_param(dep)
            #print(self._asteval.symtable)
            self._asteval.symtable[name] = par.value
            #self._asteval.symtable[name.replace('_comma_','').replace('Model', '').replace('Peaks','')] = par.value
            updated_tracker.discard(name)
            for i in updated_tracker:
                i_ = i.replace('_comma_','').replace('Model', '').replace('Peaks','')
                updated_tracker_.add(i_)

        for name in requires_update:
            _update_param(name)
            #print('classe nova; name: {}, index: {}'.format(name[:-1], int(name[-1:])))
            name_ = name.split('_comma_')
            if self.main_class:
                #print('estou auqi')
                self.main_class.upd_func_pars_ff()
                self.main_class.upd_func_pars_sf()
                self.main_class.calculate(name_[0], int(name_[2]), True, 'Model')    