## Tutorial para instalação dos pacotes necessários

Neste ponto do curso você já deve saber como abrir o terminal da sua máquina
virtual/laptop. Não assumirei que vocês sabem como usá-lo, mas adianto que 
em tempos de automação de trabalho, todos deveriam saber alguma coisa. 

- [Link com comandos básicos](http://ringo.ams.stonybrook.edu/index.php/Unix)
- [Link com tutorial de shell scripting](http://ringo.ams.stonybrook.edu/index.php/BASH_scripting)

### 1 - Instalar o Miniconda

Acesse o site `https://docs.conda.io/en/latest/miniconda.html` e faça o download 
da versão compatível com sua máquina virtual. 

Como Ubuntu é uma distribuição Linux, desça a barra de rolagem e escolha na coluna
Python 3.7 o link `Miniconda3 Linux 64-bit`. Um script `Miniconda3-py37_4.11.0-Linux-x86_64.sh` 
será salvo na pasta `Downloads`.

No terminal, digite:

``` 
cd ~/
mkdir local
mv ~/Downloads/Miniconda3-py37_4.11.0-Linux-x86_64.sh ~/local
cd local
bash Miniconda3-py37_4.11.0-Linux-x86_64.sh
```

Esses comandos te levarão ao seu diretório principal (`cd ~/`) onde você criará
uma pasta chamada `local` para onde você moverá o script salvo em `~/Downloads`.
O comando `bash Miniconda3-py37_4.11.0-Linux-x86_64.sh` instalará o Miniconda. 
Aceite todas as condições que aparecerem na tela.

**Atenção:** Toda vez que quiser usar os pacotes instalados pelo `conda` você 
deverá garantir que um ambiente conda está sendo usado. Para isso, digite:

```
which python
```
se o terminal retornar algo como `/usr/bin/python` sem referenciar `conda` ou 
`miniconda`, você deve executar o seguinte comando:

```
conda activate
```

### 2 - Instalar os pacotes usando o Miniconda.

`conda` é um administrador de pacotes em Python. Cada pacote contém livrarias de 
códigos que você pode reusar à vontade. Não é necessário reinventar a roda! Para
instalar os pacotes principais faça, na seguinte ordem:

```
conda install numpy
conda install scipy
conda install matplotlib
conda install rdkit
conda install -c conda-forge ambertools
conda install -c omnia openmm
conda install -c omnia mdtraj
```

Aceite todas as condições que aparecerem na tela em cada uma dessas instalações.

### 3 - Instalar pacotes não-incluídos no Miniconda

Acesse o site [https://github.com/m3g/packmol/releases] e faça o download da 
versão mais recente do programa `packmol` no arquivo `tar.gz`.

Para acessar o conteúdo e compilar o programa na sua máquina virtual, copie 
`packmol.tar.gz` de `~/Downloads` para `~/local` e descomprima o arquivo com o 
seguinte comando:

```
tar -xvf packmol.tar.gz
```

Esse comando criará uma pasta chamada `packmol` dentro de `~/local`. Digite:

```
cd packmol
make
```

Esse comando deve instalar `packmol` sem problemas na sua máquina virtual. 
Tendo dificuldades, confira como pode resolver o problema [neste link](http://leandro.iqm.unicamp.br/m3g/packmol/userguide.shtml#comp)

Digite no terminal

```
cd ~/local
git clone https://github.com/ParmEd/ParmEd
```
Essa sequência de comandos criará uma pasta chamada `ParmEd` em `local`. Para 
entrar e instalar o programa temos que:

```
cd ParmEd
python setup.py install
```

**Pronto! Agora você tem tudo disponível para começar a estudar a parte prática
da dinâmica molecular!**
