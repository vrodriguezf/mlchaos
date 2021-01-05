import fastai
from fastai.imports import *
from IPython.display import Audio, display, HTML, Javascript, clear_output # from tsai
import importlib

def beep(inp=1, duration=.1, n=1):
    rate = 10000
    mult = 1.6 * inp if inp else .08
    wave = np.sin(mult*np.arange(rate*duration))
    for i in range(n): 
        display(Audio(wave, rate=10000, autoplay=True))
        time.sleep(duration / .1)
        
def m_reload(package_name):
    for k,v in sys.modules.items():
        if k.startswith(package_name):
            importlib.reload(v)