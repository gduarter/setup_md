# Preparando simulações de moléculas pequenas

Temos que definir bem o tipo de simulação que será realizada:
- É uma molécula de soluto em uma caixa de solvente?
- É uma mistura de solventes?
- É uma molécula de soluto em uma mistura de solventes?

## Uma molécula de soluto em uma caixa de água

Primeiramente necessitamos a estrutura tridimensional do soluto em um arquivo `.mol2`
ou um arquivo `.pdb`. A estrutura 3D pode ser desenhada em um programa específico
(e.g., UCSF Chimera, Avogadro) ou pode ser criada a partir de SMILES, uma representação
unidimensional da estrutura da molécula, usando os pacotes OpenBabel ou RDKit.

Em posse do arquivo com a estrutura molecular, é necessário atribuir os parâmetros do
campo de forças. Isso é feito com o programa `antechamber`, do pacote `Ambertools`.
```
antechamber -i molname.pdb -fi pdb -o molname.mol2 -fo mol2 -at gaff2 -c bcc -rn LIG -nc 0
```
`-i` diz respeito ao _input_ com a estrutura da molécula de interesse; `-fi` sinaliza
o formato do _input_; `-o`, o nome do _output_; `-fo`, o formato do _output_. `-at`
determina o campo de força a ser usado; `-c`, que tipo de cargas atômicas serão
atribuídas a cada átomo. `-nc` corresponde a carga total e `-rn` é o código de 3
letras que o usuário dá a molécula. `molname` corresponde ao nome do soluto.

O programa `parmchk2` é usado para verificar se tudo está em ordem com os parâmetros
e predizê-los em casos não bem definidos.
```
 parmchk2 -i molname.mol2 -f mol2 -o molname.frcmod
```

A molécula é solvatada e os parâmetros do campo de força são aplicados no sistema
completo pelo programa `tleap`, que gerará dois arquivos `system.prmtop` e
`system.inpcrd`, contendo a topologia e as coordenadas do sistema.

O input do `tleap` é simples:

```
cat <<EOF > tl.in
source leaprc.gaff2
source leaprc.water.tip3p
LIG = loadmol2 molname.mol2
loadamberparams molname.frcmod
solvatebox LIG TIP3PBOX 12.0
setbox LIG centers 0.0
saveamberparm LIG molname.prmtop molname.inpcrd
quit
EOF

tleap -f tl.in
```
Essas linhas criam o arquivo `tl.in` e executam o programa `tleap`, que solvata
o soluto em questão com o modelo `TIP3P` de água e cria dois arquivos, `molname.prmtop`
e `molname.inpcrd`. Para criar os arquivos `.gro` e `.top` apropriados para o GROMACS,
use o script [`amber_to_gmx.py`](https://github.com/gduarter/setup_md/blob/master/zzz.scripts/amber_to_gmx.py).


## Uma molécula de soluto em um solvente diferente da água

Quando o solvente não é água, o procedimento para gerar uma caixa de simulação
exige o uso de um software específico, `packmol`, que pode ser obtido neste
[link](http://leandro.iqm.unicamp.br/m3g/packmol/download.shtml). O usuário deve
providenciar as estruturas em formato `.pdb` dos componentes da mistura e preparar
um arquivo de _input_ com os comandos abaixo:

```
cat <<EOF > input.inp
tolerance 1.5 # tolerance distance
output solute_in_solvent.pdb # output file name
filetype pdb # output file type
#
structure soluto.pdb
number NSOLU # Number of molecules
resnumbers 3 # Sequential numbering
inside cube 0. 0. 0. 30. # x, y, z coordinates of box, and length of box in Angstroms
add_amber_ter
end structure
#
structure solvente.pdb
number NSOLV # Number of molecules
resnumbers 3 # Sequential numbering
inside cube 0. 0. 0. 30.
add_amber_ter
end structure
EOF

packmol < input.inp
```

A partir dos arquivos `.pdb` de cada molécula, podemos criar os arquivos `.mol2`
de cada molécula com `antechamber` e verificar os parâmetros com o `parmchk2`,
assim como foi no caso anterior. O `tleap` deve ser rodado logo em seguida para
criar os arquivos `.prmtop` e `.inpcrd`:

```
cat <<EOF > tl.in
source leaprc.gaff2

MOL = loadmol2 solute.mol2
loadamberparams solute.frcmod
SOL = loadmol2 solvent.mol2
loadamberparams solvent.frcmod

fullbox = loadPdB solute_in_solvent.pdb

setbox fullbox centers

saveAmberParm fullbox solute_in_solvent.prmtop solute_in_solvent.inpcrd
quit
EOF

tleap -f tl.in
```

Os arquivos `.inpcrd` e `.prmtop` podem ser transformados em arquivos para o
GROMACS da mesma forma que foi mencionado na sessão anterior.
