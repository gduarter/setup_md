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
e `molname.inpcrd`. Para criar os arquivo `.gro` e `.top` apropriados para o GROMACS,
use o script `amber_to_gmx.py`.


## Uma molécula de soluto em um solvente diferente da água
