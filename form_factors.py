import numpy as np
import lmfit as lmf
import scipy.special as sc 

##### Variables:

#Essa variável contém os modelos de fatores de forma implemetados
all = ['Background', 'Sphere', 'Core Shell Sphere', 'Spheroid', 'Core Shell Spheroid', 'Cylinder', 'Core Shell Cylinder']

##### Functions:
def ruido(obj, bg=0, a=1):
    '''Esta função simula um ruído e um backgroud no fator de forma do objeto'''
    if a != 0:
        i_ = (obj.intensidade+bg)/a
    obj.intensidade = i_+(np.random.randn(len(i_))/np.sqrt(i_))
    return

##### Classes:

class FormFactor:
    '''Essa é uma classe genérica de fatores de forma. Ela contém as funções que são iguais para todos os modelos. A demais classes de fatores de forma devem herdar os atributos dessa classe'''
    def __init__(self):
        '''Essa é a função de inicialização da classe. São declaradas algumas variáveis que podem ser necessárias em outras funções.
        A variável "self.dist" é utulizada para guardar qual distribuição é aplicada.
        A variável "self.Rg" guarda o valor calculado do raio de giro.
        A varável "self.current_method" guarda qual método escolhido para o input dos parâmetros do modelo.'''
        self.dist = None
        self.Rg = None #raio de giro guinier
        #self.real_value = None
        self.current_method = 0

    def setq(self, q): 
        '''Esta função atribui o valor de q para o cálculo do ff. A variável "q" deve ser um array.'''
        self.q = q #vetor de espalhamento q
    
    def ff3 (self, q, func_pars): 
        '''Essa função calcula o fator de forma usando os valores recebidos para cada parâmetro.
        A variável "q" deve ser um array com os pontos onde o fator de forma deve ser calculado.
        A variável "func_pars" deve ser um dicionário dos parâmetros do modelo.
        A variável "self.intensidade" será um array com a mesma dismensão que "q" com os calculados da intensidade de espalhamento para cada ponto.'''
        self.update_all_params(q, func_pars)
        self.intensidade = self.ff()
        return self.intensidade

    def update_all_params(self, q, func_pars):
        '''Esta função atuliza os valores dos parâmetros salvos no objeto.
        A variável "q" deve ser um array com os pontos onde o fator de forma deve ser calculado.
        A variável "func_pars" deve ser um dicionário dos parâmetros do modelo.'''
        for i in func_pars:
            exec('self.'+i+r"= func_pars[i]['value']")
        self.q = q

    def ff_dist(self, sigma, funcpars, p, d, min_amp = -3, max_amp = 3, n = 50):
        '''Esta função retorna o fator de forma considerando uma distribuição de valores em um dos parâmetros.
        O valor médio da distribuição será o valor na salvo da variavel no objeto
        O parâmetro 'sigma' é o valor sigma da distribuição em percentual
        A variável "d" definie qual distribuição deve ser aplicada.
        A variável "func_pars" deve ser um dicionário dos parâmetros do modelo.
        O parâmetro 'p' é sobre qual variável a distribuição se aplicada
        A o parâmetro 'amp' diz para até quanto vezes o valor de sigma distante da média a distribuição seré calculada
        '''
        x, y = self.distribuicao(sigma, funcpars, p, d, min_amp, max_amp, n) # x e y da distribuição.
        norm = (np.sum(y)) #normalização
        intensidade = []
        if len(self.methods)>1 and p>0:
            for i in x: #cálculo do fator de forma para cada valor da distribuição aplicada  
                i_ = self.func(i, self.current_method, False) 
                #print('i_: ', i_)
                intensidade.append(self.ff2(i_, funcpars, p)) #calculo da intensidade
        else:
            for i in x: #cálculo do fator de forma para cada valro da distribuição aplicada  
                intensidade.append(self.ff2(i, funcpars, p)) #calculo da intensidade
        intensidade = np.transpose(intensidade)
        a_ = (intensidade*y)/norm #integrando
        i_dist = np.sum(a_,axis=1) #integração
        return i_dist
    
    def distribuicao(self, sigma, funcpars, p, d, min_amp = -3, max_amp = 3, n = 50): 
        '''Esta função calcula a distribuição de tamanhos de uma determinada variável.
        O valor médio da distribuição será o valor na salvo da variavel no objeto
        O parâmetro 'sigma' é o valor sigma da distribuição em percentual
        O parâmetro 'p' é sobre qual variável a distribuição se aplicada
        A o parâmetro 'amp' diz para até quanto vezes o valor de sigma distante da média a distribuição seré calculada'''
        from distribution import Distribuicao
        self.param_index = p
        param = self.dist_param(funcpars, p)
        self.dist = Distribuicao() #calcula a distribuição
        #print('68: ', sigma, param, d, min_amp, max_amp, n)
        x, y = self.dist.distribuicao(sigma, param, d, min_amp, max_amp, n)#calcula a distribuição
        return x, y

'''
Algumas variáveis serão necessárias em todas as classes de fatores de forma. Elas estarão listadas abaixo:
    - A variável "self.func_pars" é o dicionário com os valores dos paramêtros do modelo. Cada item do dicionário corresponde a um parâmetro. Cada item também deve ser um dicionário, 
    contendo os itens:
        -'value': valor do parâmetro (float);
        -'min': valor mínimo do parâmetro (float);
        -'máx': valor máximo do parâmetro (float);
        -'symbol': Texto exibido para identificar o parâmetro na GUI. (string).
        
    - A variável "self.methods" é um array com os possíveis diferentes parâmetros de input dos modelos.

    - A variável "self.changeparam" com os parâmetros que serão trocados quando o método de input é mudado. 
        essa é uma possível mudanda a ser feita. Da maneira que está agora apenas um slider pode ser modificado quando o método de input é trocado.

    - A variável "self.dist_params" é um discionário com os parâmetros do modelo em que é possível aplicar a distribuição de tamanhos.
'''

class Background(FormFactor):
    '''Essa classe é utilizada para adiconar um background aos dados. Ela foi escrita como uma classe igual as demais classes de fatores de forma.'''
    def __init__(self, bg = 0):
        '''A variável "bg" deve ser do tipo float. O valor dela será o valor do background adicionado.'''
        super().__init__()
        self.name = 'Background'
        self.bg = bg

        self.func_pars = {'bg':{'value': bg, 'min': 0, 'max': 500, 'symbol': 'bg'}}

        self.methods = ['Bg']
        self.changeparam = ['bg']
        self.dist_params = {0:['bg']}

    def ff(self):
        '''Essa função retonar um array com a mesma dismensão da variável do vetor de espalhamento "self.q" com o valor definido do background na variável "self.bg".'''
        self.intensidade = np.zeros(len(self.q))+self.bg
        return self.intensidade
    
    def ff_classic(self, q, bg):
        '''Essa função retonar um array com a mesma dismensão da variável do vetor de espalhamento "self.q" com o valor definido do background na variável "self.bg".
        A variável "q" deve ser um array com os valores de interesse do vetor de espalhamento;
        A variável "bg" deve ser do tipo float com o valor desejado do background.'''
        intensidade = np.zeros(len(q))+bg
        return intensidade

class Esfera(FormFactor):
    ''' Esta classe calcula o fator de forma de uma esfera'''
    def __init__(self, r = 50, rho = 1):
        '''Essa função inicializa a classe "Esfera".
        A variável "r" deve ser um float com o valor desejado do raio da esfera.
        A variável "rho" deve ser um float com o valor desejado da diferença entre a densidade eletronica do meio e da esfera.'''
        super().__init__()
        self.name = 'Esfera'
        self.rho = rho #densidade eletronica
        self.r = r #raio

        #Dicionário com os valores máximos, mínimos e atuais dos sliders dos parâmetros necessários para o cálculo do fator de forma da esfera
        self.func_pars = {'r':{'value':r, 'min':1, 'max':100, 'symbol':'R'},
                          'rho':{'value':rho, 'min':-2, 'max':2, 'symbol':'\u03C1'}}
        
        self.methods = ['R'] #Array com as diferentes maneiras aceitas para o input dos parâmetros
        self.changeparam = ['r'] #Array com as keys do dicionários de sliders que devarão ser modificadas ao trocar a variável de input
        self.dist_params = {0:['R']} #Dicionário com os parâmetros em que é possível aplicar uma distribuição

    def rg_t(self):
        '''Esta função retorna o raio de giro teórico da esfera'''
        return np.sqrt(3/5)*self.r 
    
    def ff(self): 
        '''Esta função calcula o fator de forma da esfera com os valores dos parâmetros que já estão salvos no objeto'''
        self.intensidade = (((4*(np.pi)*self.rho* ((np.sin(self.q*self.r)-self.q*self.r*np.cos(self.q*self.r)) / (self.q**3)))**2))/((4/3)*np.pi*self.r**3) #fator de forma da esfera
        return self.intensidade
    
    def ff_classic(self, q, r, rho):
        '''Esta função retorna o fator de forma de uma esfera com raio r e densidade eletronica rho dados.
        A variável "q" deve ser um array com os valores de interesse do vetor de espalhamento;
        A variável "r" deve ser um float com o valor desejado do raio da esfera.
        A variável "rho" deve ser um float com o valor desejado da diferença entre a densidade eletronica do meio e da esfera.'''
        intensidade = ((4*(np.pi)*rho* ((np.sin(q*r)-q*r*np.cos(q*r)) / (q**3)))**2)/((4/3)*np.pi*r**3)
        return intensidade
    
    def ff2(self, r, funcpars, param_index): 
        '''Esta função retorna o fator de forma da esfera para um dado raio. Para os outros parâmetros usa os valores já salvos no objetos.
        Ela é utulizada para calcular o fator de forma considerando uma distribuição de tamanhos aplicada no raio.
        A variável "r" deve ser um float com o valor desejado do raio da esfera.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "param_index" indica em qual parâmetro a distribuição será aplicada. Nesse caso ela só poder ser aplicada no raio, então essa variável não é utilizada.
        Ela foi mantida nessa função pois é necessária em no caso de outros fatores de forma.'''
        rho = funcpars['rho']['value']
        return ((4*(np.pi)*rho*((np.sin(self.q*r)-self.q*r*np.cos(self.q*r))/(self.q**3)))**2) / ((4/3)*np.pi*r**3) #fator de forma da esfera

    def dist_param(self, funcpars, p):
        '''Essa função retorna o valor do parâmetro que se deseja aplicar a distribuição.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "p" indica em qual parâmetro a distribuição será aplicada. Nesse caso ela só poder ser aplicada no raio, então essa variável não é utilizada.
        Ela foi mantida nessa função pois é necessária em no caso de outros fatores de forma.
        '''
        self.p = 'r'
        return funcpars['r']['value'] #valor médio
    
class Elipsoide(FormFactor):
    '''Esta classe calcula o fator de forma de um elipsoide de revolução'''
    def __init__(self, r = 50, nu = 1.2, rho = 1):
        '''Essa função inicializa a classe "Elipsoide".
        A variável "r" deve ser um float com o valor de um dos eixos da elipse que quando rotacionada gera o elipsoide.
        A variável "nu" deve ser um float com o valor desejado para a razão entre os eixos da elipse que quando rotacionada gera o elipsoide.
        A variável "rho" deve ser um float com o valor desejado da diferença entre a densidade eletronica do meio e do elipsoide.'''
        super().__init__()
        self.name = 'Elipsóide'
        self.rho = rho #densidade eletronica
        self.nu = nu #nu = R2/R1
        self.r = r #raio 1

        #Dicionário com os valores máximos, mínimos e atuais dos sliders dos parâmetros necessário para o cálculo do fator de forma de um elipsóide de revolução
        self.func_pars = {'r':{'value':r, 'min':1, 'max':100, 'symbol':'R'}, 'nu':{'value':nu, 'min':0, 'max':5, 'symbol':'\u03BD'}, 
                          'rho':{'value':rho, 'min':-2, 'max':2, 'symbol':'\u03C1'}}
        
        self.methods = ['R, \u03BD', 'R\u2081, R\u2082', 'R, t+R'] #Array com as diferentes maneiras aceitas para o input dos parâmetros

        #Array com as key do dicionário referente ao slider devará ser modificado ao trocar a variavel de input e o novo texto que deverá ser exibido:
        self.changeparam = ['nu', '\u03BD'] 

        #Dicionário com os parâmetros em que é possível aplicar uma distribuição. Os números são referente aos diferentes métodos de input
        self.dist_params = {0:['R', '\u03BD'], 1:['R\u2081', 'R\u2082'], 2:['R', 't']} 

        #array com os angulos para a integração, de forma que a intensidade calculada será a media do objeto orientado em todas as direções    
        self.theta = np.linspace(0, np.pi, 400)

    def rg_t(self):
        '''Esta função retorna o valor do raio de giro teórico do elipsóide'''
        return np.sqrt((2*self.r**2+(self.nu*self.r)**2)/5)
    
    def ff(self):
        '''Esta função calcula o fator de forma de um elipsóide de revolução com os valores dos parâmetros que já estão salvos na objeto'''
        self.theta = np.linspace(0, np.pi, 181) #valores de theta
        raio_ = self.r*np.sqrt((self.nu**2)*(np.cos(self.theta))**2 + (np.sin(self.theta))**2)
        norm = sum(np.sin(self.theta)) #normalização para a integral numerica
        Int = []
        for i in range(len(self.theta)): #cálculo do fator de forma para cada alpha
            Int.append((3* ((np.sin(self.q*raio_[i])-self.q*raio_[i]*np.cos(self.q*raio_[i])) / ((raio_[i]*self.q)**3)))**2) #cálculo dos pontos para integração
        a = np.sin(self.theta)*np.transpose(Int)/norm #integrando
        volume = (4./3.)*self.nu*np.pi*self.r**3
        self.intensidade = (np.sum(a, axis=1)*(volume*self.rho)**2)/volume #integração
        return self.intensidade
    
    def ff_classic(self, q, r1, r2, rho):
        '''Esta função calcula o fator de forma de um elipsóide de revolução com os valores dos parâmetros informados.
        A variável "q" deve ser um array com os valores de interesse do vetor de espalhamento;
        A variável "r1" deve ser um float com o valor de um dos eixos da elipse que quando rotacionada gera o elipsoide.
        A variável "n2" deve ser um float com o valor outro eixo da elipse que quando rotacionada gera o elipsoide.
        A variável "rho" deve ser um float com o valor desejado da diferença entre a densidade eletronica do meio e do elipsoide.'''
        theta = np.linspace(0, np.pi, 181) #valores de theta
        r = r1
        nu = r2/r1
        raio_ = r*np.sqrt((nu**2)*(np.cos(theta))**2 + (np.sin(theta))**2)
        norm = sum(np.sin(theta)) #normalização para a integral numerica
        Int = []
        for i in range(len(theta)): #cálculo do fator de forma para cada alpha
            Int.append((3* ((np.sin(q*raio_[i])-q*raio_[i]*np.cos(q*raio_[i])) / ((raio_[i]*q)**3)))**2) #cálculo dos pontos para integração
        a = np.sin(theta)*np.transpose(Int)/norm #integrando
        volume = (4./3.)*nu*np.pi*r**3
        intensidade = (np.sum(a, axis=1)*(volume*rho)**2)/volume #integração
        return intensidade
    
    def ff2(self, par, funcpars, param_index):
        '''Esta função é utilizada para o calculo do fator de forma quando uma distribuição de tamanhos é aplicada.
        Com ela é possível calcular o fator de forma fornecendo o valor de um dos parâmetros e usar para os ooutros os valores já salvos no objeto.
        A variável "par" deve ser do tipo float e será o valor do parâmetro em que a distribuição está sendo aplicada.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "param_index" indica em qual parâmetro a distribuição será aplicada.'''
        rho = funcpars['rho']['value']
        nu = funcpars['nu']['value']
        if param_index == 0:
            r = par
            nu = funcpars['nu']['value']
        else:#if param_index == 1:
            r = funcpars['r']['value']
            nu = par

        raio_ = r*np.sqrt((nu**2)*(np.cos(self.theta))**2 + (np.sin(self.theta))**2)
        norm = sum(np.sin(self.theta))
        Int = []
        for i in range(len(self.theta)):
            Int.append((3* ((np.sin(self.q*raio_[i])-self.q*raio_[i]*np.cos(self.q*raio_[i])) / ((raio_[i]*self.q)**3)))**2)
        a = np.sin(self.theta)*np.transpose(Int)/norm
        volume = (4./3.)*nu*np.pi*r**3
        intensidade = (np.sum(a, axis=1)*(volume*rho)**2)/volume
        return intensidade
    
    def func(self, a , i, atribute = True): 
        '''Esta função transforma os inputs R2 e t+R em nu, para o calculo do fator de forma.
        A variável 'a' é o input
        A variável 'i' informa qual método está sendo utilizado
        A variável "atribute" informa se a função deve trocar o metodo de input salvo na variável "self.current_method". Ela existe para que seja possível utilizar essa função para
        fazer a conversão dos valores entre os diferentes métodos sem alterar os dados salvos no objeto.'''
        if atribute:
            self.current_method = i
        if i == 0: #neste caso input ja é o valor de nu
            self.changeparam = ['nu', '\u03BD']
            return a #o valor retornado é o próprio input
        if i == 1: #neste caso o input é R2, assim nu = R2/R1
            self.changeparam = ['nu', 'R\u2082']
            return a/self.slider['r'].value() #nu = R2/R1
        if i == 2: #neste caso o input é t. R2 = R1+t -> nu = (R1+t)/R1
            self.changeparam = ['nu', 't']
            return (a+self.slider['r'].value())/self.slider['r'].value() #nu = (R1+t)/R1

    def value_methods(self, method, slider):
        '''Esta função transforma o valor de nu em t e R2 , para o cálculo do fator de forma do elipsóde de revolução.
        A variavel 'slider' é um dicionário contendo os sliders da user interface.
        A variavel 'method' informa qual método esta sendo utilizado.
        A variavel 'self.value' é um dicionário que guarda os valores mínimos, máximos e atuais dos parâmetros.
        Essa função é utilizada para calcular os valores que devem ser exibidos no slider da GUI quando o método de input é trocado.'''
        self.slider = slider #salva os sliders como uma variavel da classe
        if method == 0: #nu
            self.value = {'value': self.func_pars['nu']['value'], 'min':self.func_pars['nu']['min'], 'max':self.func_pars['nu']['max']}
        if method == 1:  #R2 = nu*R1
            value = self.func_pars['nu']['value']*self.func_pars['r']['value']
            min = self.func_pars['nu']['min']*self.func_pars['r']['value']
            max = self.func_pars['nu']['max']*self.func_pars['r']['value']
            self.value = {'value':value, 'min':min, 'max':max}
        if method == 2:  #t+R1 = nu*R1 -> t = nu*R1 - R1
            value = self.func_pars['nu']['value']*self.func_pars['r']['value']-self.func_pars['r']['value']
            min = self.func_pars['nu']['min']*self.func_pars['r']['min']-self.func_pars['r']['min']
            max = self.func_pars['nu']['max']*self.func_pars['r']['max']-self.func_pars['r']['max']
            self.value = {'value':value, 'min':min, 'max':max}
    
    def dist_param(self, funcpars, p):
        '''Essa função retorna o valor do parâmetro que se deseja aplicar a distribuição.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "p" indica em qual parâmetro a distribuição será aplicada.'''
        if p == 0:
            self.p = 'r'
            return funcpars['r']['value']
        else:
            self.p = 'nu'
            return self.slider['nu'].value()

class Cilindro(FormFactor):
    '''Esta classe calcula o fator de forma de um cilindro'''
    def __init__(self, r = 50, l = 120, rho = 1):
        '''Essa função inicializa a classe "Cilindro".
        A variável "r" deve ser um float com o valor desejado do raio do cilindro.
        A variável "l" deve ser um float com o valor desejado do comprimento do clilindro.
        A variável "rho" deve ser um float com o valor desejado da diferença entre a densidade eletronica do meio e do cilindro'''
        super().__init__()
        self.name = 'Cilindro'
        self.rho = rho #densidade eletronica
        self.r = r #raio
        self.l = l #comprimento

        #Dicionário com os valores máximos, mínimos e atuais dos sliders dos parâmetros necessário para o cálculo do fator de forma de um cilindro
        self.func_pars = {'r':{'value':r, 'min':1, 'max':100, 'symbol':'R'}, 'l':{'value':l, 'min':1, 'max':300, 'symbol':'L'}, 
                          'rho':{'value':rho, 'min':-2, 'max':2, 'symbol':'\u03C1'}}
        
        self.methods = ['R, L'] #Array com as diferentes maneiras aceitas para o input dos parâmetros

        #Array com as chaves do dicionário referente ao slider devará ser modificado ao trocar a variavel de input e o novo texto que deverá ser exibido:
        self.changeparam = ['r'] 

        #Dicionário com os parâmetros em que é possível aplicar uma distribuição. Os números são referente aos diferentes métodos de input:
        self.dist_params = {0:['R', 'L']}

    def rg_t(self):
        '''Esta função retorna o valor do raio de giro teórico do Cilindro'''
        return np.sqrt(((self.r**2/2)+(self.l**2/12)))

    def ff(self):
        '''Esta função calcula o fator de forma de um cilindro com os valores dos parâmetros que já estão salvos na objeto'''
        intensity = []
        self.alpha = np.linspace(0, np.pi/2, 181) #valores de alpha para a integração
        volume = np.pi*self.r**2*self.l  
        for i in range(len(self.alpha)): #cálculo do fator de forma para cada alpha
            x = 0.5*self.q*self.l*np.cos(self.alpha[i])
            y = self.q*self.r*np.sin(self.alpha[i])
            if np.cos(self.alpha[i]) == 0:
                intensity.append((sc.jv(1,y)/y)**2)
            elif np.sin(self.alpha[i])==0:
                intensity.append((0.5*np.sin(x)/x)**2)
            else:
                intensity.append((np.sin(x)/x*sc.jv(1,y)/y)**2) #calculo da intensidade para cada alpha
        norm = sum(np.sin(self.alpha)) #normalização da integração
        a = (2*self.rho*volume)**2*np.sin(self.alpha)*np.transpose(intensity)/norm #integrando
        self.intensidade = np.sum(a, axis=1)/volume #calculo da integral
        return self.intensidade
    
    def ff_classic(self, q, r, l, rho):
        '''Esta função retorna o fator de forma de um cilindro com raio r, comprimento l e densidade eletronica rho dados.
        A variável "q" deve ser um array com os valores de interesse do vetor de espalhamento;
        A variável "r" deve ser um float com o valor desejado do raio do cilindro.
        A variável "l" deve ser um float com o valor desejado do comprimento do cilindro.
        A variável "rho" deve ser um float com o valor desejado da diferença entre a densidade eletronica do meio e do cilindro.'''
        intensity = []
        alpha = np.linspace(0, np.pi/2, 180) #valores de alpha para a integração
        volume = np.pi*r**2*l  
        for i in range(len(alpha)): #cálculo do fator de forma para cada alpha
            x = 0.5*q*l*np.cos(alpha[i])
            y = q*r*np.sin(alpha[i])
            if np.cos(alpha[i]) == 0:
                intensity.append((sc.jv(1,y)/y)**2)
            elif np.sin(alpha[i])==0:
                intensity.append((0.5*np.sin(x)/x)**2)
            else:
                intensity.append((np.sin(x)/x*sc.jv(1,y)/y)**2) #calculo da intensidade para cada alpha
        norm = sum(np.sin(alpha)) #normalização da integração
        a = (2*rho*volume)**2*np.sin(alpha)*np.transpose(intensity)/norm #integrando
        intensidade = np.sum(a, axis=1)/volume #calculo da integral
        return intensidade
    
    def ff2(self, par, funcpars, param_index):#, r=0,l=0): 
        '''Esta função é utilizada para o calculo do fator de forma quando uma distribuição de tamanhos é aplicada.
        Com ela é possível calcular o fator de forma fornecendo o valor de um dos parâmetros e usar para os outros os valores já salvos no objeto.
        A variável "par" deve ser do tipo float e será o valor do parâmetro em que a distribuição está sendo aplicada.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "param_index" indica em qual parâmetro a distribuição será aplicada.'''
        intensity = []
        rho = funcpars['rho']['value']
        if param_index==0:
            r=par
            l = funcpars['l']['value']
        else:#if param_index == 1:
            r = funcpars['r']['value']
            l = par
        volume = np.pi*r**2*l  
        for i in range(len(self.alpha)):
            x = 0.5*self.q*l*np.cos(self.alpha[i])
            y = self.q*r*np.sin(self.alpha[i])
            if np.cos(self.alpha[i]) == 0:
                intensity.append((sc.jv(1,y)/y)**2)
            elif np.sin(self.alpha[i])==0:
                intensity.append((0.5*np.sin(x)/x)**2)
            else:
                intensity.append((np.sin(x)/x*sc.jv(1,y)/y)**2) 
        norm = sum(np.sin(self.alpha))
        a = (2*rho*volume)**2*np.sin(self.alpha)*np.transpose(intensity)/norm
        intensidade = np.sum(a, axis=1)/volume
        return intensidade

    def dist_param(self, funcpars, p):
        '''Essa função retorna o valor do parâmetro que se deseja aplicar a distribuição.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "p" indica em qual parâmetro a distribuição será aplicada.'''
        if p == 0:
            self.p = 'r'
            return funcpars['r']['value']
        else:
            self.p = 'l'
            return funcpars['l']['value']
            
class CascaEsferica(FormFactor):
    '''Esta classe calcula o fator de forma de uma Casca Esférica'''
    def __init__(self, rc = 50, rs = 60, rho1 = 1, rho2 = 1):
        '''Essa função inicializa a classe "Esfera".
        A variável "rc" deve ser um float com o valor do raio interno da casca esférica.
        A variável "ns" deve ser um float com o valor do raio externo da casca esférica.
        A variável "rho1" deve ser um float com o valor desejado da diferença entre a densidade eletronica do núcleo e da casca.
        A variável "rho2" deve ser um float com o valor desejado da diferença entre a densidade eletronica do casca e do meio.'''
        super().__init__()
        self.name = 'Casca Esférica'
        self.rho1 = rho1 #rho_core - rho_shell
        self.rho2 = rho2 #rho_shell - rho_sol
        self.rc = rc #raio_core
        self.rs = rs #raio_shell

        #Dicionário com os valores máximos, mínimos e atuais dos sliders dos parâmetros necessário para o cálculo do fator de forma de uma Casca Esférica:
        self.func_pars = {'rc':{'value':rc, 'min':1, 'max':100, 'symbol':'R\u2081'}, 'rs':{'value':rs, 'min':1, 'max':200, 'symbol':'R\u2082'},
                           'rho1':{'value':rho1, 'min':-2, 'max':2, 'symbol':'\u03C1\u2081'}, 'rho2':{'value':rho2, 'min':-2, 'max':2,
                            'symbol':'\u03C1\u2082'}} 
        
        self.methods = ['R\u2081, R\u2082', 'R, \u03BDR', 'R, t+R'] #Array com as diferentes maneiras aceitas para o input dos parâmetros

        #Array com as keys do dicionário referente ao slider devará ser modificado ao trocar a variavel de input e o novo texto que deverá ser exibido:
        self.changeparam = ['rs', 'R\u2082']

        #Dicionário com os parâmetros em que é possível aplicar uma distribuição. Os números são referente aos diferentes métodos de input:
        self.dist_params = {0:['R\u2081', 'R\u2082'], 1:['R', '\u03BD'], 2:['R', 't']} 

    def rg_t(self):
        '''Esta função retorna o valor do raio de giro teórico de uma casca esférica'''
        return np.sqrt((3/5)*(self.rs**5-self.rc**5)/(self.rs**3-self.rc**3))
   
    def func(self, a , i, atribute = True):
        '''Esta função transforma os inputs 'nu' e 't' em 'R2', para o cálculo do fator de forma da casca esferérica
        A variável 'a' o input
        A variável 'i' informa qual método está sendo utilizado.
        A variável "atribute" informa se a função deve trocar o metodo de input salvo na variável "self.current_method". Ela existe para que seja possível utilizar essa função para
        fazer a conversão dos valores entre os diferentes métodos sem alterar os dados salvos no objeto.'''
        if atribute:  
            self.real_value = a          
            self.current_method = i
        if i == 0: #neste caso input ja é o valor de R2
            self.changeparam = ['rs', 'R\u2082']
            return a
        if i == 1:  #neste caso o input é nu, assim R2 = nu*R1
            self.changeparam = ['rs', '\u03BD']
            return a*self.slider['rc'].value() #R2 = nu*R1
        if i == 2: #neste caso o input é t. R2 = R1+t
            self.changeparam = ['rs', 't']
            return a+self.slider['rc'].value() #R2 = R1+t

    def value_methods(self, method, slider):
        '''Esta função transforma o valor de R2 em nu e t para o cálculo do fator de forma de uma casca esférica
        A variavel 'slider' é um dicionário contendo os sliders da user interface
        A variavel 'method' informa qual método esta sendo utilizado
        A variavel 'self.value' é um dicionário guarda os valores minimos, máximos e atuais dos parâmetros.
        Essa função é utilizada para calcular os valores que devem ser exibidos no slider da GUI quando o método de input é trocado.'''
        self.slider = slider #salva os sliders como uma variavel da classe
        if method == 0: #R2
            self.value = {'value': self.func_pars['rs']['value'], 'min':self.func_pars['rs']['min'], 'max':self.func_pars['rs']['max']}
        if method == 1: #R2 = nu*R1 -> nu = R2/R1
            self.value = {'value': self.func_pars['rs']['value']/self.func_pars['rc']['value'],
                          'min':self.func_pars['rs']['min']/self.func_pars['rc']['value'],
                          'max':self.func_pars['rs']['max']/self.func_pars['rc']['value']}
        if method == 2: #R2 = t+R1 -> t = R2-R1
            self.value = {'value': self.func_pars['rs']['value']-self.func_pars['rc']['value'],
                          'min':self.func_pars['rs']['min']-self.func_pars['rc']['value'],
                          'max':self.func_pars['rs']['max']-self.func_pars['rc']['value']}

    def ff(self):  
        '''Esta função calcula o fator de forma de uma casca esférica com os valores dos parâmetros que já estão salvos na objeto.'''
        vs = (4/3)*np.pi*(self.rs**3) #Volume shell
        vc = (4/3)*np.pi*self.rc**3 #volume core
        ks = (np.sin(self.q*self.rs) - self.q*self.rs*np.cos(self.q*self.rs))/(self.q*self.rs)**3
        kc = (np.sin(self.q*self.rc) - self.q*self.rc*np.cos(self.q*self.rc))/(self.q*self.rc)**3
        self.intensidade = (((3)*(vc*self.rho1*kc - vs*self.rho2*ks))**2)/vs
        return self.intensidade

    def ff_classic(self, q, rc, rs, rho1, rho2):
        '''Esta função calcula o fator de forma de uma casca esférica com os valores dos parâmetros informados.
        A variável "rc" deve ser um float com o valor do raio interno da casca esférica.
        A variável "ns" deve ser um float com o valor do raio externo da casca esférica.
        A variável "rho1" deve ser um float com o valor desejado da diferença entre a densidade eletronica do núcleo e da casca.
        A variável "rho2" deve ser um float com o valor desejado da diferença entre a densidade eletronica do casca e do meio.'''
        vs = (4/3)*np.pi*(rs**3) #Volume shell
        vc = (4/3)*np.pi*rc**3 #volume core
        ks = (np.sin(q*rs) - q*rs*np.cos(q*rs))/(q*rs)**3
        kc = (np.sin(q*rc) - q*rc*np.cos(q*rc))/(q*rc)**3
        intensidade = (((3)*(vc*rho1*kc - vs*rho2*ks))**2)/vs
        return intensidade
    
    def ff2(self, r, funcpars, param_index): 
        '''Esta função é utilizada para o calculo do fator de forma quando uma distribuição de tamanhos é aplicada.
        Com ela é possível calcular o fator de forma fornecendo o valor de um dos parâmetros e usar para os ooutros os valores já salvos no objeto.
        A variável "par" deve ser do tipo float e será o valor do parâmetro em que a distribuição está sendo aplicada.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "param_index" indica em qual parâmetro a distribuição será aplicada.'''
        rho1 = funcpars['rho1']['value']
        rho2 = funcpars['rho2']['value']
        if param_index == 0:
            rc = r
            rs = funcpars['rs']['value']

        else:#if param_index == 1:
            rs = r
            rc = funcpars['rc']['value']

        vs = (4/3)*np.pi*(rs**3)
        vc = (4/3)*np.pi*rc**3
        ks = (np.sin(self.q*rs) - self.q*rs*np.cos(self.q*rs))/(self.q*rs)**3
        kc = (np.sin(self.q*rc) - self.q*rc*np.cos(self.q*rc))/(self.q*rc)**3
        intensidade = (((3)*(vc*rho1*kc - vs*rho2*ks))**2)/vs
        return intensidade
    
    def dist_param(self, funcpars, p):
        '''Essa função retorna o valor do parâmetro que se deseja aplicar a distribuição.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "p" indica em qual parâmetro a distribuição será aplicada.'''
        #self.ff3(self.q, self.func_pars)
        if p == 0:
            self.p = 'rc'
            return funcpars['rc']['value']
        else:
            self.p = 'rs'  
            return self.slider['rs'].value()

class ElipsoideCoreShell(FormFactor):
    '''Esta classe calcula o fator de forma de um elipsóide de revolução core shell.'''
    def __init__(self, r_e = 20, x_c = 3, t_s = 30, x_s=1, rho1 = 0.5, rho2 = 0.5):
        '''Essa função inicializa a classe "ElipsoideCoreShell".
        A variável "r_e" deve ser um float com o valor do raio equatorial interno do elipsoide .
        A variável "x_c" deve ser um float com o valor desejado para a razão entre o raio equatorial interno e o raio polar interno.
        A variável "t_s" deve ser um float com o valor da espessura da casca no equador.
        A variável "x_s" deve ser um float com o valor da razaõ entre a espessura da casca no polo e no equador.
        A variável "rho1" deve ser um float com o valor desejado da diferença entre a densidade eletronica do núcleo e da casca.
        A variável "rho2" deve ser um float com o valor desejado da diferença entre a densidade eletronica do casca e do meio.'''
        super().__init__()
        self.name = 'Elipsóide Core Shell'
        self.r_e, self.x_c, self.t_s, self.x_s = r_e, x_c, t_s, x_s   #raio equatorial, axial ratio of core, X = r_polar/r_equatorial, thickness of shell at equator, ratio of thickness of shell at pole to that at equator
        self.rho1, self.rho2 = rho1, rho2   #rho1 = rho_core-rho_shell, rho2 = rho_shell-rho_solv

        #Dicionário com os valores máximos, mínimos e atuais dos sliders dos parâmetros necessário para o cálculo do fator de forma de um elipsóde core shell
        self.func_pars = {'r_e':{'value':r_e, 'min':1, 'max':100, 'symbol':'R\u2091'}, 'x_c':{'value':x_c, 'min':0.1, 'max':10, 'symbol':'X_c'}, 
                          't_s':{'value':t_s, 'min':1e-3, 'max':50, 'symbol':'t\u209b'}, 'x_s':{'value':x_s, 'min':0.1, 'max':10, 'symbol':'X\u209b'},
                          'rho1':{'value':rho1, 'min':-2, 'max':2, 'symbol':'\u03C1\u2081'}, 'rho2':{'value':rho2, 'min':-2, 'max':2, 'symbol':'\u03C1\u2082'}}
        
        self.methods = ['R\u209A, X_c, t\u209b, '] #Array com as diferentes maneiras aceitas para o input dos parâmetros

        #Array com as chaves do dicionário referente ao slider devará ser modificado ao trocar a variavel de input e o novo texto que deverá ser exibido:
        self.changeparam = ['r_p'] 

        #Dicionário com os parâmetros em que é possível aplicar uma distribuição. Os números são referente aos diferentes métodos de input:
        self.dist_params = {0:['R\u209a', 'X_c', 't_s', 'X\u209b']}
    
    def rg_t(self):
        '''Esta função retorna o valor do raio de giro teórico do Elpsoide cor shell'''
        return None
    
    def ff(self):
        '''Esta função calcula o fator de forma de um elipsóide core shell com os valores dos parâmetros que já estão salvos na objeto'''
        intensity = []
        self.mu = np.linspace(0, np.pi/2, 91) #valores de mu para a integração

        def j1(x):
            return (np.sin(x)-x*np.cos(x))/x**3
        
        def _x(r_e, r_p, alpha):
            return np.sqrt(r_e**2*np.sin(alpha)**2+r_p**2*np.cos(alpha)**2)


        for mu in self.mu: #cálculo do fator de forma para cada mu
            #print('mu: ', mu)
            v_c = (4/3)*np.pi*self.r_e**3*self.x_c
            v_s = (4/3)*np.pi*(self.r_e + self.t_s)**2*(self.r_e*self.x_c+self.t_s*self.x_s)

            _x_c = self.q*_x(self.r_e, self.r_e*self.x_c, mu)
            #print('x_c:', x_c)
            _x_s = self.q*_x(self.r_e+self.t_s, self.r_e*self.x_c +self.t_s*self.x_s, mu)
            #print('x_t:', x_t)
            intensity.append((self.rho1*v_c*3*j1(_x_c) + self.rho2*v_s*3*j1(_x_s))**2) #calculo da intensidade para cada alpha
        norm = sum(np.sin(self.mu)) #normalização da integração
        a = np.transpose(intensity)*np.sin(self.mu)/norm #integrando
        self.intensidade = np.sum(a, axis=1)/v_s  #calculo da integral
        return self.intensidade
    
    def ff_classic(self, q, r_e = 20, x_c = 3, t_s = 30, x_s=1, rho1 = 0.5, rho2 = 0.5):
        '''Esta função calcula o fator de forma de um elipsóide de revolução core shell com os valores dos parâmetros informados.
        A variável "r_e" deve ser um float com o valor do raio equatorial interno do elipsoide .
        A variável "x_c" deve ser um float com o valor desejado para a razão entre o raio equatorial interno e o raio polar interno.
        A variável "t_s" deve ser um float com o valor da espessura da casca no equador.
        A variável "x_s" deve ser um float com o valor da razaõ entre a espessura da casca no polo e no equador.
        A variável "rho1" deve ser um float com o valor desejado da diferença entre a densidade eletronica do núcleo e da casca.
        A variável "rho2" deve ser um float com o valor desejado da diferença entre a densidade eletronica do casca e do meio.
        '''
        intensity = []
        _mu = np.linspace(0, np.pi/2, 91) #valores de mu para a integração

        def j1(x):
            return (np.sin(x)-x*np.cos(x))/x**3
        
        def _x(r_e, r_p, alpha):
            return np.sqrt(r_e**2*np.sin(alpha)**2+r_p**2*np.cos(alpha)**2)


        for mu in _mu: #cálculo do fator de forma para cada mu
            v_c = (4/3)*np.pi*r_e**3*x_c #volume core
            v_s = (4/3)*np.pi*(r_e + t_s)**2*(r_e*x_c+t_s*x_s) #volume shell

            _x_c = q*_x(r_e, r_e*x_c, mu)
            _x_s = q*_x(r_e+t_s, r_e*x_c +t_s*x_s, mu)
            intensity.append((rho1*v_c*3*j1(_x_c) + rho2*v_s*3*j1(_x_s))**2) #calculo da intensidade para cada alpha
        norm = sum(np.sin(_mu)) #normalização da integração
        a = np.transpose(intensity)*np.sin(_mu)/norm #integrando
        intensidade = np.sum(a, axis=1)/v_s #calculo da integral
        return intensidade

    def ff2(self, par, funcpars, param_index):
        '''Esta função é utilizada para o calculo do fator de forma quando uma distribuição de tamanhos é aplicada.
        Com ela é possível calcular o fator de forma fornecendo o valor de um dos parâmetros e usar para os ooutros os valores já salvos no objeto.
        A variável "par" deve ser do tipo float e será o valor do parâmetro em que a distribuição está sendo aplicada.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "param_index" indica em qual parâmetro a distribuição será aplicada.'''
        intensity = []
        rho1 = funcpars['rho1']['value']
        rho2 = funcpars['rho2']['value']

        #verificação de qual parâmetro o valor está sendo dado
        if param_index==0:
            r_e=par
            x_c = funcpars['x_c']['value']
            t_s = funcpars['t_s']['value']
            x_s = funcpars['x_s']['value']      
        elif param_index == 1:
            r_e = funcpars['r_e']['value']
            x_c = par
            t_s = funcpars['t_s']['value']
            x_s = funcpars['x_s']['value']
        elif param_index == 2:
            r_e = funcpars['r_e']['value']
            x_c = funcpars['x_c']['value']
            t_s = par
            x_s = funcpars['x_s']['value']
        else:
            r_e = funcpars['r_e']['value']
            x_c = funcpars['x_c']['value']
            t_s = funcpars['t_s']['value']
            x_s = par
        
        intensity = []
        def j1(x):
            return (np.sin(x)-x*np.cos(x))/x**3
        
        def _x(r_e, r_p, alpha):
            return np.sqrt(r_e**2*np.sin(alpha)**2+r_p**2*np.cos(alpha)**2)
        for mu in self.mu: #cálculo do fator de forma para cada mu
            v_c = (4/3)*np.pi*r_e**3*x_c
            v_s = (4/3)*np.pi*(r_e + t_s)**2*(r_e*x_c+t_s*x_s)

            _x_c = self.q*_x(r_e, r_e*x_c, mu)
            _x_s = self.q*_x(r_e+t_s, r_e*x_c +t_s*x_s, mu)
            intensity.append((rho1*v_c*3*j1(_x_c) + rho2*v_s*3*j1(_x_s))**2) #calculo da intensidade para cada alpha
        norm = sum(np.sin(self.mu)) #normalização da integração
        a = np.transpose(intensity)*np.sin(self.mu)/norm #integrando
        intensidade = np.sum(a, axis=1)/v_s #calculo da integral
        return intensidade
    
    def dist_param(self, funcpars, p):
        '''Essa função retorna o valor do parâmetro que se deseja aplicar a distribuição.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "p" indica em qual parâmetro a distribuição será aplicada.'''
        if p == 0:
            self.p = 'r_e'
            return funcpars['r_e']['value']
        elif p == 1:
            self.p = 'x_c'
            return funcpars['x_c']['value']
        elif p == 2:
            self.p = 't_s'
            return funcpars['t_s']['value']
        else:
            self.p = 'x_s'
            return funcpars['x_s']['value']

class CilindroCoreShell(FormFactor):
    '''Esta classes calcula o fator de forma de um Cilindro Core Shell'''
    def __init__(self, r = 20, t1 = 20, t2=20, l = 400, rho1 = 0, rho2 = 3):
        '''Essa função inicializa a classe "CilindroCoreShell".
        A variável "r" deve ser um float com o valor do raio interno do cilindro.
        A variável "t1" deve ser um float com o valor da diferença entre o raio externo e interno do cilindro.
        a variável "t2" deve ser um float com o valor da diferença entre o comprimento interno e externo do cilindro.
        a variável "l" deve ser um float com o valor do comprimento interno so cilindro.
        A variável "rho1" deve ser um float com o valor desejado da diferença entre a densidade eletronica do núcleo e da casca.
        A variável "rho2" deve ser um float com o valor desejado da diferença entre a densidade eletronica do casca e do meio.'''
        super().__init__()
        self.name = 'Cilindro Core Shell'
        self.r = r #raio do cilindro
        self.t1 = t1  #espessura da casca
        self.t2 = t2 
        self.l = l #comprimento
        self.rho1 = rho1 #rho_core-rho_shell
        self.rho2 = rho2 #rho_sell-rho_solv

        #Dicionário com os valores máximos, mínimos e atuais dos sliders dos parâmetros necessário para o cálculo do fator de forma de um Cilindro Core SHell:
        self.func_pars = {'r':{'value':r, 'min':1, 'max':100, 'symbol':'R'}, 't1':{'value':t1, 'min':1, 'max':100, 'symbol':'t1'}, 't2':{'value':t2, 'min':1, 'max':100, 'symbol':'t2'},
                           'l':{'value':l, 'min':1, 'max':500, 'symbol':'l'}, 'rho1':{'value':rho1, 'min':-5, 'max':5, 'symbol':'\u03C1\u2081'},
                           'rho2':{'value':rho2, 'min':-5, 'max':5,'symbol':'\u03C1\u2082'}} 
        
        self.methods = ['R, t\u2081, t\u2082, l'] #Array com as diferentes maneiras aceitas para o input dos parâmetros

        #Array com as key do dicionário referente ao slider devará ser modificado ao trocar a variavel de input e o novo texto que deverá ser exibido:
        self.changeparam = ['r']

        #Dicionário com os parâmetros em que é possível aplicar uma distribuição. Os números são referente aos diferentes métodos de input:
        self.dist_params = {0:['R', 't\u2081', 't\u2082', 'l']} 

    def rg_t(self):
        '''Esta função retorna o valor do raio de giro teórico do Cilindro Core Shell'''
        return None
    
    def ff(self):
        '''Esta função calcula o fator de forma de um Cilindro Core Shell com os valores dos parâmetros que já estão salvos na objeto'''
        intensity = []
        self.alpha = np.linspace(0, np.pi/2, 181) #valores de alpha para a integração
        volume_core = np.pi*self.r**2*self.l
        volume_shell = np.pi*(self.r+self.t1)**2*(self.l+2*self.t2)  
        for i in range(len(self.alpha)): #cálculo do fator de forma para cada alpha
            x_c = 0.5*self.q*self.l*np.cos(self.alpha[i])
            y_c = self.q*self.r*np.sin(self.alpha[i])
            x_s = self.q*(0.5*self.l+self.t2)*np.cos(self.alpha[i])
            y_s = self.q*(self.r+self.t1)*np.sin(self.alpha[i])
            if np.cos(self.alpha[i]) == 0:
                intensity.append(((2*self.rho1*volume_core*sc.jv(1,y_c)/y_c)+(2*self.rho2*volume_shell*sc.jv(1,y_s)/y_s))**2)
            elif np.sin(self.alpha[i])==0:
                intensity.append(((2*self.rho1*volume_core*0.5*np.sin(x_c)/x_c)+(2*self.rho2*volume_shell*0.5*np.sin(x_s)/x_s))**2)
            else:
                intensity.append(((2*self.rho1*volume_core*np.sin(x_c)/x_c*sc.jv(1,y_c)/y_c)+(2*self.rho2*volume_shell*np.sin(x_s)/x_s*sc.jv(1,y_s)/y_s))**2) #calculo da intensidade para cada alpha
        norm = sum(np.sin(self.alpha)) #normalização da integração
        a = np.sin(self.alpha)*np.transpose(intensity)/norm #integrando
        self.intensidade = np.sum(a, axis=1)/volume_shell #calculo da integral
        return self.intensidade

    def ff_classic(self, q, r, t1, t2, l, rho1, rho2):
        '''Esta função calcula o fator de forma de um cilindro core sehll com os valores dos parâmetros informados.
        A variável "r" deve ser um float com o valor do raio interno do cilindro.
        A variável "t1" deve ser um float com o valor da diferença entre o raio externo e interno do cilindro.
        a variável "t2" deve ser um float com o valor da diferença entre o comprimento interno e externo do cilindro.
        a variável "l" deve ser um float com o valor do comprimento interno so cilindro.
        A variável "rho1" deve ser um float com o valor desejado da diferença entre a densidade eletronica do núcleo e da casca.
        A variável "rho2" deve ser um float com o valor desejado da diferença entre a densidade eletronica do casca e do meio.'''
        intensity = []
        alpha = np.linspace(0, np.pi/2, 181) #valores de alpha para a integração
        volume_core = np.pi*r**2*l
        volume_shell = np.pi*(r+t1)**2*(l+2*t2)  
        for i in range(len(alpha)): #cálculo do fator de forma para cada alpha
            x_c = 0.5*q*l*np.cos(alpha[i])
            y_c = q*r*np.sin(alpha[i])
            x_s = q*(0.5*l+t2)*np.cos(alpha[i])
            y_s = q*(r+t1)*np.sin(alpha[i])
            if np.cos(alpha[i]) == 0:
                intensity.append(((2*rho1*volume_core*sc.jv(1,y_c)/y_c)+(2*rho2*volume_shell*sc.jv(1,y_s)/y_s))**2)
            elif np.sin(alpha[i])==0:
                intensity.append(((2*rho1*volume_core*0.5*np.sin(x_c)/x_c)+(2*rho2*volume_shell*0.5*np.sin(x_s)/x_s))**2)
            else:
                intensity.append(((2*rho1*volume_core*np.sin(x_c)/x_c*sc.jv(1,y_c)/y_c)+(2*rho2*volume_shell*np.sin(x_s)/x_s*sc.jv(1,y_s)/y_s))**2) #calculo da intensidade para cada alpha
        norm = sum(np.sin(alpha)) #normalização da integração
        a = np.sin(alpha)*np.transpose(intensity)/norm #integrando
        intensidade = np.sum(a, axis=1)/volume_shell #calculo da integral
        return intensidade

    def ff2(self, par, funcpars, param_index):
        '''Esta função é utilizada para o calculo do fator de forma quando uma distribuição de tamanhos é aplicada.
        Com ela é possível calcular o fator de forma fornecendo o valor de um dos parâmetros e usar para os ooutros os valores já salvos no objeto.
        A variável "par" deve ser do tipo float e será o valor do parâmetro em que a distribuição está sendo aplicada.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "param_index" indica em qual parâmetro a distribuição será aplicada.'''
        intensity = []
        rho1 = funcpars['rho1']['value']
        rho2 = funcpars['rho2']['value']

        #verificação de qual parâmetro o valor está sendo dado
        if param_index==0:
            r=par
            t1 = funcpars['t1']['value']
            t2 = funcpars['t2']['value']
            l = funcpars['l']['value']
        elif param_index == 1:
            r=funcpars['r']['value']
            t1 = par
            t2 = funcpars['t2']['value']
            l = funcpars['l']['value']
        elif param_index == 2:
            r=funcpars['r']['value']
            t1 = funcpars['t1']['value']
            t2 = par
            l = funcpars['l']['value']
        else:
            r=funcpars['r']['value']
            t1 = funcpars['t1']['value']
            t2 = funcpars['t2']['value']
            l = par

        volume_core = np.pi*r**2*l
        volume_shell = np.pi*(r+t1)**2*(l+2*t2)  
        for i in range(len(self.alpha)): #cálculo do fator de forma para cada alpha
            x_c = 0.5*self.q*l*np.cos(self.alpha[i])
            y_c = self.q*r*np.sin(self.alpha[i])
            x_s = self.q*(0.5*l+t2)*np.cos(self.alpha[i])
            y_s = self.q*(r+t1)*np.sin(self.alpha[i])
            if np.cos(self.alpha[i]) == 0:
                intensity.append(((2*rho1*volume_core*sc.jv(1,y_c)/y_c)+(2*rho2*volume_shell*sc.jv(1,y_s)/y_s))**2)
            elif np.sin(self.alpha[i])==0:
                intensity.append(((2*rho1*volume_core*0.5*np.sin(x_c)/x_c)+(2*rho2*volume_shell*0.5*np.sin(x_s)/x_s))**2)
            else:
                intensity.append(((2*rho1*volume_core*np.sin(x_c)/x_c*sc.jv(1,y_c)/y_c)+(2*rho2*volume_shell*np.sin(x_s)/x_s*sc.jv(1,y_s)/y_s))**2) #calculo da intensidade para cada alpha
        norm = sum(np.sin(self.alpha)) #normalização da integração
        a = np.sin(self.alpha)*np.transpose(intensity)/norm #integrando
        intensidade = np.sum(a, axis=1)/volume_shell #calculo da integral
        return intensidade

    def dist_param(self, funcpars, p):
        '''Essa função retorna o valor do parâmetro que se deseja aplicar a distribuição.
        A variável "funpars" é o dicionário com os valores dos parâmetros do modelo.
        A variável "p" indica em qual parâmetro a distribuição será aplicada.'''
        if p == 0:
            self.p = 'r'
            return funcpars['r']['value']
        elif p == 1:
            self.p = 't1'
            return funcpars['t1']['value']
        elif p == 2:
            self.p = 't2'
            return funcpars['t2']['value']
        elif p == 3:
            self.p = 'l'
            return funcpars['l']['value']

class MultiCascaEsferica:
    "Essa classe calcula o fator de forma de uma esfera com diversas camada. Essa classe não possui a funções necessárias para o uso na GUI."
    def __init__(self, r, rho, q):
        self.rho = rho #vetor com as diferencas das densidades eletronicas
        self.r = r  #vetor com os valores de raio de cada camada     
        self.q = q
        self.Rg = None
        self.Rg_t = None

        self.ff()
    def ff(self):
        v, k = [], []
        for i in range(len(self.r)):
            #v.append((4/3)*np.pi*(self.r[i]**3))
            k.append((4/3)*np.pi*(self.r[i]**3)*self.rho[i]*(np.sin(self.q*self.r[i]) - self.q*self.r[i]*np.cos(self.q*self.r[i]))/(self.q*self.r[i])**3)
        self.intensidade = (3*np.sum(k, axis=0))**2
        return self.intensidade

class RaioDeGiro: 
    '''Esta classe calcula o raio de giro pelo método de guinier de um determinado objeto'''
    def __init__(self, obj): #recebe uma classe fator de forma
        obj.qRgmax = 1 #q*Raio_giro máximo
        self.obj = obj #objeto fator de forma
        # variável 'raio' = array com j raios de giro calculádos por guinier
        #'r_sq' = array com os valores R^2 dos ajustes lineares para cada valor de raio de giro calculado
        #'k' = indice do maior valor de q usado
        #'j' = número de raios de giro calculados
        #'coef_l' = coeficiente linear do ajuste para linear
        #'coef_a' = coeficiente angular do ajuste
        raio, r_sq, k, j, coef_l, coef_a = self.raio_de_giro(obj.q, obj.intensidade) 
        obj.Rgs = raio #salva os raios de giro calculados na objeto
        obj.R2 = r_sq #salva os valores de R^2 no objeto
        obj.index_qmax = k #salva o índice k no objeto
        obj.index_ = j #salva o indice j no objeto
        obj.coef_l = coef_l #salva o coeficiente linear no objeto
        obj.coef_a = coef_a #salva o coeficiente angular no objeto
        obj.q_max = obj.q[k-1] #salva o valor máximo de q no objeto
        obj.Rg = raio[j-1] #salva o valor do raio de giro no objeto
        return None
    
    def setff(self, ff):
        '''Esta função atribui um fator de forma ao objeto e calcula seu raio de giro'''
        self.obj = ff
        self.__init__(self.obj)

    def set_qRg(self, valor): 
        '''Esta função muda o valor máximo de q*Rg. Por padrão, q*Rg = 1'''
        self.qRgmax = valor
        self.__init__(self.obj)

    def raio_de_giro(self, q1, Int): 
        '''Esta função calcula o raio de giro por guinier '''
        #calcula o primeiro valor do raio de giro para k=10
        k = 10 
        i2=Int[:k] #intensidade ao quadrado
        q2=q1[:k]**2   #q ao quadrado
        log_i = np.log(i2) #logaritimo da intensidade ao quadrado

        Rgs=[]
        R2 = []
        
        gmod = lmf.models.ExpressionModel("coef_a*x + coef_l")
        result = gmod.fit(log_i, x=q2, coef_a = 1, coef_l = 0)

        Rgs.append(np.sqrt(-3*result.params['coef_a'].value))  #coef_a = -(Rg^2)/3

        j=0
        #Calculo do raio de giro ate qRg
        while q1[k]*Rgs[j]<self.obj.qRgmax and k<len(q1)-1:
            k = k+1
            j+=1
            i2=Int[:k]   #intensidade ao quadrado
            q2=q1[:k]**2   #q ao quadrado

            log_i = np.log(i2) #logaritimo da intensidade

            result = gmod.fit(log_i, x=q2, coef_a = 1, coef_l = 0)
            Rgs.append(np.sqrt(-3*result.params['coef_a'].value))  #coef_a = -(Rg^2)/3
            
            #R2.append(model.score(q2.reshape((-1, 1)), log_i)) #Coeficiente R^2 do ajuste
        return Rgs, R2, k, j, result.params['coef_l'].value, result.params['coef_a'].value
    
class FileData: 
    '''Esta classe cria um objeto utilizando um arquivo. '''
    def __init__(self, x, y, yerr = None, xerr = None, q_unit = None):
        self.intensidade = y #intensidade importada
        self.q = x #vetor espalhamento do arquivo
        self.q_original = x
        self.yerr = yerr
        self.xerr = xerr
        self.Rg = None #raio de giro por guinier
        self.path = None
        self.q_unit = q_unit

dic = {0:Background, 1:Esfera, 2:CascaEsferica, 3:Elipsoide, 4:ElipsoideCoreShell, 5:Cilindro, 6:CilindroCoreShell}