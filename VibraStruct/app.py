"""Interactive browser dashboard. Launch: python -m streamlit run app.py."""
import io
import json
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from vibrastruct import sdof, receptance, chain, modes, rayleigh, response, energy, harmonic_response

st.set_page_config(page_title='VibraStruct | Vibration Analysis',page_icon='〰️',layout='wide')
st.title('VibraStruct')
st.caption('STRUCTURAL DYNAMICS & VIBRATION ANALYSIS')
st.write('Explore resonance, damping and coupled motion. Change the inputs to compare the response of idealised mechanical systems.')


def csv_download(label, values, header, name):
    buf=io.StringIO(); np.savetxt(buf,values,delimiter=',',header=header,comments='')
    st.download_button(label,buf.getvalue(),file_name=name,mime='text/csv')


def figure_download(fig,name):
    buf=io.BytesIO(); fig.savefig(buf,format='png',dpi=180,bbox_inches='tight')
    st.download_button('Download figure',buf.getvalue(),file_name=name,mime='image/png')


with st.sidebar:
    st.header('Analysis settings')
    model=st.radio('System',['Single mass','Coupled three-mass system'])
    st.caption('Inputs use SI units. All defaults are synthetic examples.')
    with st.form('parameters'):
        if model=='Single mass':
            mass=st.number_input('Mass [kg]',10.,10000.,250.,10.)
            stiffness=st.number_input('Stiffness [N/m]',10000.,100000000.,1000000.,10000.)
        else:
            masses=[]; springs=[]
            for i,(mv,kv) in enumerate(zip([250.,180.,120.],[1000000.,800000.,600000.])):
                masses.append(st.number_input(f'Mass {i+1} [kg]',10.,10000.,mv,10.))
                springs.append(st.number_input(f'Spring {i+1} [N/m]',10000.,100000000.,kv,10000.))
        zeta=st.slider('Damping ratio [%]',0.5,20.,3.,0.5)/100
        amplitude=st.number_input('Harmonic force amplitude [N]',0.,10000.,100.,10.)
        duration=st.slider('Free-vibration duration [s]',0.5,10.,3.,0.5)
        displacement=st.number_input('Initial displacement [mm]',-20.,20.,2.,0.5)
        if model!='Single mass':
            force_station=st.selectbox('Force applied at mass',[1,2,3],index=2)
            initial_station=st.selectbox('Initially displaced mass',[1,2,3],index=2)
        st.form_submit_button('Update analysis',type='primary')
    st.caption('Free vibration starts from rest. Harmonic forcing is a separate steady-state study.')

if model=='Single mass':
    wn,c=sdof(mass,stiffness,zeta); M=np.array([[mass]]); K=np.array([[stiffness]]); C=np.array([[c]])
    w=np.array([wn]); P=np.array([[1/np.sqrt(mass)]]); F=np.array([amplitude]); x0=np.array([displacement/1000])
    modal_z=np.array([zeta])
else:
    M,K=chain(masses,springs); w,P=modes(M,K); C=rayleigh(M,K,zeta)
    F=np.zeros(3); F[force_station-1]=amplitude
    x0=np.zeros(3); x0[initial_station-1]=displacement/1000
    modal_z=np.diag(P.T@C@P)/(2*w)

freqs=w/(2*np.pi); n=len(w)
cols=st.columns(n+1)
for j,a in enumerate(freqs): cols[j].metric(f'Natural frequency {j+1}',f'{a:.3f} Hz')
cols[-1].metric('Damping target',f'{100*zeta:.1f}%')
frequency_tab,time_tab,mode_tab,method_tab=st.tabs(['Frequency response','Free vibration','Mode shapes','Methods & assumptions'])
with frequency_tab:
    st.subheader('Harmonic displacement response')
    st.caption('Automatic frequency range includes all natural frequencies. Dotted lines mark the undamped natural frequencies.')
    f=np.linspace(0,max(30.,1.6*freqs[-1]),2401)
    H=harmonic_response(M,K,C,2*np.pi*f,F)
    fig,ax=plt.subplots(figsize=(10,4))
    for j in range(n): ax.plot(f,np.abs(H[:,j])*1000,label=f'Mass {j+1}')
    for a in freqs: ax.axvline(a,color='grey',ls=':',alpha=0.5)
    ax.set(xlabel='Excitation frequency [Hz]',ylabel='Displacement amplitude [mm]'); ax.grid(alpha=0.2); ax.legend()
    st.pyplot(fig); figure_download(fig,'vibrastruct_frequency.png'); plt.close(fig)
    csv_download('Download frequency response CSV',np.c_[f,np.abs(H),np.angle(H)],'frequency_Hz,'+','.join(f'amplitude{j+1}_m' for j in range(n))+','+','.join(f'phase{j+1}_rad' for j in range(n)),'frequency_response.csv')
with time_tab:
    st.subheader('Motion after release')
    st.caption('Only the initial displacement acts in this study; there is no external force after release.')
    # Output and internal steps resolve the fastest natural period with at least 40 samples.
    dt=min(0.001,1/(40*freqs[-1])); count=int(np.ceil(duration/dt))+1
    t=np.linspace(0,duration,count)
    @st.cache_data(show_spinner='Calculating the transient response…',max_entries=8)
    def transient(M,K,C,t,x0):
        return response(M,K,C,t,x0,np.zeros(len(x0)))
    x,v=transient(M,K,C,t,x0); E=energy(M,K,x,v)
    stride=max(1,len(t)//6000)
    fig,axs=plt.subplots(2,1,figsize=(10,6),sharex=True)
    for j in range(n): axs[0].plot(t[::stride],x[::stride,j]*1000,label=f'Mass {j+1}')
    axs[0].set_ylabel('Displacement [mm]'); axs[0].legend(); axs[0].grid(alpha=0.2)
    axs[1].plot(t[::stride],E[::stride],color='#139c88'); axs[1].set(xlabel='Time [s]',ylabel='Mechanical energy [J]'); axs[1].grid(alpha=0.2)
    fig.tight_layout(); st.pyplot(fig); figure_download(fig,'vibrastruct_transient.png'); plt.close(fig)
    csv_download('Download transient CSV',np.c_[t,x,v,E],'time_s,'+','.join(f'x{j+1}_m' for j in range(n))+','+','.join(f'v{j+1}_m_s' for j in range(n))+',energy_J','transient.csv')
with mode_tab:
    st.subheader('Relative motion of each mass')
    fig,axs=plt.subplots(1,n,figsize=(10,4),squeeze=False,sharey=True)
    for j,ax in enumerate(axs[0]):
        shape=P[:,j]/np.max(np.abs(P[:,j])); ax.plot(np.r_[0,shape],np.arange(n+1),'o-')
        ax.set(xlabel='Relative displacement',title=f'Mode {j+1} · {freqs[j]:.2f} Hz',yticks=np.arange(n+1)); ax.grid(alpha=0.2)
    axs[0,0].set_ylabel('Mass station (0 = ground)'); fig.tight_layout()
    st.pyplot(fig); figure_download(fig,'vibrastruct_modes.png'); plt.close(fig)
    st.caption('Shapes are scaled to unit maximum amplitude. Connecting lines indicate mass stations, not a continuous beam deformation.')
    st.dataframe({'Mode':np.arange(1,n+1),'Frequency [Hz]':freqs,'Modal damping [%]':modal_z*100},hide_index=True)
    csv_download('Download modal results CSV',np.c_[np.arange(1,n+1),freqs,modal_z,P.T],'mode,frequency_Hz,damping_ratio,'+','.join(f'phi_mass{j+1}' for j in range(n)),'modes.csv')
with method_tab:
    st.markdown('''### Linear structural dynamics
The model solves **M ẍ + C ẋ + K x = f(t)**. Natural modes solve the generalized eigenvalue problem **K φ = ω² M φ**. Frequencies are reported in Hz. Harmonic amplitudes use a cosine forcing convention.

The coupled chain is grounded at spring 1 and free beyond mass 3. Each spring connects consecutive stations. Small translations, constant masses, linear springs and viscous damping are assumed. Rayleigh damping matches the selected ratio at the lowest and highest modes; the middle ratio can differ.

Transient motion uses adaptive RK45 integration. Results are sampled at least 40 times per fastest natural period. The frequency grid provides sampled peaks; very lightly damped peaks may need a finer sweep for precise amplitudes.

### Scope
Inputs are synthetic examples. The model is educational and has no experimental calibration. It does not calculate stress, model a continuum structure or determine safe operating speeds for real machinery.

### Verification
The accompanying Python package includes analytical benchmarks, energy checks and independent coupled-system comparisons. Run `python -m unittest discover -s tests -v` to check the installed package.''')
    payload={'model':model,'mass_matrix_kg':M.tolist(),'stiffness_matrix_N_m':K.tolist(),'damping_matrix_Ns_m':C.tolist(),'force_amplitude_N':F.tolist(),'initial_displacement_m':x0.tolist(),'duration_s':duration,'natural_frequencies_Hz':freqs.tolist(),'modal_damping_ratios':modal_z.tolist()}
    st.download_button('Download model settings JSON',json.dumps(payload,indent=2),file_name='vibrastruct_settings.json',mime='application/json')
st.caption('VibraStruct · SI units · Linear models · Reproducible numerical analysis')
