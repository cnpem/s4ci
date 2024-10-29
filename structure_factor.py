import numpy as np

##############################################
##           FATORES DE ESTRUTURA           ##
##############################################
##### Variables:
all = ['Hard Sphere']
aprox = ['DA', 'LMA']#, 'PSF']

##### Classes:
class StructureFacture:
    '''Essa é uma classe genérica de fatores de estrutura. Ela contém as funções que são iguais para todos os modelos. 
    A demais classes de fatores de estrutura devem herdar os atributos dessa classe.'''
    def setq(self, q): 
        '''Esta função atribui o valor de q para o cálculo do sf'''
        self.q = q #vetor de espalhamento q

    def update_all_params(self, q, func_pars):
        '''Esta função atuliza os valores dos parâmetros salvos no objeto.
        A variável "q" deve ser um array com os pontos onde o fator de estrutura deve ser calculado.
        A variável "func_pars" deve ser um dicionário dos parâmetros do modelo.'''
        for i in func_pars:
            exec('self.'+i+r"= func_pars[i]['value']")
        self.q = q

    def sf3(self, q, func_pars):
        '''Essa função calcula o fator de forma usando os valores recebidos para cada parâmetro.
        A variável "q" deve ser um array com os pontos onde o fator de forma deve ser calculado.
        A variável "func_pars" deve ser um dicionário dos parâmetros do modelo.
        A variável "self.intensidade" será um array com a mesma dismensão que "q" com os valores calculados do fator de estrutura para cada ponto.'''
        self.update_all_params(q, func_pars)
        self.intensidade = self.sf()
        return self.intensidade

    def set_Ref(self, r_ef):
        '''Essa classe atribui um valor do raio efetivo para o modelo.
        #TODO Eu acredito que essa função na verdade não deva ficar aqui. Preciso ver.'''
        self.r_ef = r_ef

class HardSphere(StructureFacture):
    '''Esta classe calcula o fator de estrutura considerando o modelo de esfera rígida.'''
    def __init__(self, r_ef, fp=0.2):
        super().__init__()
        self.r_ef = r_ef
        self.fp = fp

        self.func_pars = {'fp':{'value':fp, 'min':0, 'max':0.999, 'symbol':'fp'}, 'r_ef':{'value':r_ef, 'min':1, 'max':300, 'symbol':'R_ef'}}

        self.methods = ['R_ef'] #Array com as diferentes maneiras aceitas para o input dos parâmetros
        self.changeparam = ['r_ef'] #Array com as keys do dicionários de sliders que devarão ser modificadas ao trocar a variável de input

    def sf(self):
        '''Esta função calcula o fator de estrutura de uma esfera dura com os valores dos parâmetros que já estão salvos no objeto'''

        alpha = (1+2*self.fp)**2/(1-self.fp)**4
        beta = -6*self.fp*(1+self.fp/2)**2/(1-self.fp)**4
        gamma = self.fp*alpha/2

        a = 2*self.r_ef*self.q

        par1 = (np.sin(a) - a*np.cos(a))/a**2
        par2 = (2*a*np.sin(a) + (2-a**2)*np.cos(a) - 2)/a**3
        par3 = (-a**4*np.cos(a) + 4*((3*a**2-6)*np.cos(a) + (a**3-6*a)*np.sin(a) +6))/a**5

        g = alpha*par1 + beta*par2 + gamma*par3

        self.intensidade = 1/(1 + 24*self.fp*g/(a))

        return self.intensidade
    
    def sf2(self, r_ef):
        alpha = (1+2*self.fp)**2/(1-self.fp)**4
        beta = -6*self.fp*(1+self.fp/2)**2/(1-self.fp)**4
        gamma = self.fp*alpha/2

        a = 2*r_ef*self.q

        par1 = (np.sin(a) - a*np.cos(a))/a**2
        par2 = (2*a*np.sin(a) + (2-a**2)*np.cos(a) - 2)/a**3
        par3 = (-a**4*np.cos(a) + 4*((3*a**2-6)*np.cos(a) + (a**3-6*a)*np.sin(a) +6))/a**5

        g = alpha*par1 + beta*par2 + gamma*par3

        self.intensidade = 1/(1 + 24*self.fp*g/(a))

        return self.intensidade
    
    def sf_funcpars(self, q, funcpars):
        fp = funcpars['fp']['value']
        r_ef = funcpars['r_ef']['value']
        alpha = (1+2*fp)**2/(1-fp)**4
        beta = -6*fp*(1+fp/2)**2/(1-fp)**4
        gamma = fp*alpha/2

        a = 2*r_ef*q

        par1 = (np.sin(a) - a*np.cos(a))/a**2
        par2 = (2*a*np.sin(a) + (2-a**2)*np.cos(a) - 2)/a**3
        par3 = (-a**4*np.cos(a) + 4*((3*a**2-6)*np.cos(a) + (a**3-6*a)*np.sin(a) +6))/a**5

        g = alpha*par1 + beta*par2 + gamma*par3

        intensidade = 1/(1 + 24*fp*g/(a))

        return intensidade