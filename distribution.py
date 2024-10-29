import numpy as np
import scipy.special as sc 
from scipy import stats
import matplotlib.pyplot as plt

class Distribuicao: 
    '''Esta classe calcula uma distribuição
    A variável 'sigma' é o valor sigma sigma em percentual da distribuição
    A variável 'mu' é a média da distribuição
    A variável 'd' qual distribuição será aplicada
    A variável 'min_amp' diz qual o valor inicial para o cálculo da distribuição: x_inicial = mu + min_amp*sigma
    A variável 'max_amp' diz qual o valor final para o cálculo da distribuição: x_final = mu + max_amp*sigma
    A varipavel 'n' diz quantos pontos da distribuição serão calculados'''
    def __init__(self, sigma=0.1, mu=0, d='ln', min_amp=-3, max_amp=3, n=50):
        self.update_all(sigma, mu, d, min_amp, max_amp, n)

    def update_all(self, sigma, mu, d, min_amp, max_amp, n):
        '''Esta função atualiza os valores dos parêmtros dentro das variáveis'''
        #print('1:',sigma)
        #print('2:',mu)
        self.sigma = sigma
        self.mu = mu
        self.d = d
        self.min_amp = min_amp
        self.max_amp = max_amp
        self.n = n

    def distribuicao(self, sigma, mu, d, min_amp=-3, max_amp=3, n = 50):
        '''Esta função calcula uma distribuição
        A variável 'sigma' é o valor sigma em percentual da distribuição
        A variável 'mu' é a média da distribuição
        A variável 'd' qual distribuição será aplicada
        A variável 'min_amp' diz qual o valor inicial para o cálculo da distribuição: x_inicial = mu + min_amp*sigma
        A variável 'max_amp' diz qual o valor final para o cálculo da distribuição: x_final = mu + max_amp*sigma
        A varipavel 'n' diz quantos pontos da distribuição serão calculados'''
        self.update_all(sigma, mu, d, min_amp, max_amp, n)
        if d == 'n': #gaussiana
            sigma = sigma*mu
            self.x = np.linspace(mu+(min_amp*sigma), mu+(max_amp*sigma), n)
            self.desv_param = (1./(sigma*np.sqrt(2*np.pi))) * np.exp(-(0.5)*((self.x-mu)/sigma)**2)
            #fig, ax2 = plt.subplots()
            #ax2.plot(obj.x, obj.desv_param)        
            return self.x, self.desv_param
        
        if d == 'ln': #log-normal         
            
            #### Cálculo dos valores iniciais e finais de x da distribuição:

            # Para a log normal é encontrado os valores de x_i e x_f de forma que a distribuição seja calculada dentro do mesmo intervalo de probabilidade da gaussiana
            # A função scipy.special.ndtr retorna a probabilidade P(x<=X) da distribuição normal padrão (mu=0, sigma=1).
            # Assim obtemos os valores 'min_p' e 'max_p'. Esses valores nos dizem sob qual intervalo de probabilidades cumulativa iremos calcular a distribuição lognormal   
            min_p = sc.ndtr(((min_amp))) 
            max_p = sc.ndtr(((max_amp-sigma)))

            #A função scipy.lognorm.ppf é o inverso da função de distribuição acumulada, ou seja, para determinada probabilidade p ela retorna o valor de X tal que P(x<=X)=p.
            # Aplicando essa função aos valores de probabilidade 'min_p' e 'max_p' calculados logo acima temos o intevalo de x na qual a distribuição lognormal será calculado
            logn = stats.lognorm(sigma, scale=mu)
            min_x = logn.ppf(min_p) #valor mínimo de x
            max_x = logn.ppf(max_p) #valor máximo de x
            #################################################################

            self.x = np.linspace((min_x),(max_x), n) #valores de x

            self.desv_param = (1./(self.x*sigma*np.sqrt(2*np.pi))) * np.exp((-0.5)*((np.log(self.x/mu))/sigma)**2) #calculo da distribuição
            return self.x, self.desv_param
    
    def plot(self):
        '''Esta função plota a distribuição calculada'''
        fig, ax_plot = plt.subplots()
        ax_plot.plot(self.x, self.desv_param)
        fig.show()

    def return_dist_curve(self):
        return self.x, self.desv_param