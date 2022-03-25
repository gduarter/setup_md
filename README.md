# Tutorial de Dinâmica Molecular

Este é um tutorial para quem deseja aprender a preparar simulações de dinâmica
molecular para o software GROMACS. Como muitos alunos não têm experiência com
`bash` ou `python` _scripting_, este tutorial tem um mecanismo automatizado para
 a geração de arquivos de _input_ do GROMACS (`.gro`, `.top`, `.mdp`).

Para poder usar o material deste repositório você deve instalar alguns pacotes
em seu computador. Veja o [tutorial de instalação](zzz.markdown/install.md) dos
pacotes para tirar suas dúvidas.

## Gerando arquivos de _input_

O procedimento para gerar os arquivos é simples. Basta digitar no terminal o
seguinte comando:

```
bash organize.sh -u solutes.csv -s solvents.csv
```

Cada um dos arquivos `CSV` (_comma separated values_) contém uma tabela em que
cada molécula é identificada pela sua representação unidimensional (`SMILES`) e
pelo seu nome. Os scripts ignoram linhas começadas com o caractere `#`.

Exemplo de arquivo `CSV`:
```
SMILES,Name
c1ccccc1-c2ccccc2,biphenyl
O=C(C(C)NC(C)(C)C)c1cccc(Cl)c1,bupropion
NCCc1cc(O)c(O)cc1,dopamine
C[C@@](O)(CCO)CC(=O)O,mevalonicacid
O=C(O)[C@@H](N)Cc2c1cc(O)ccc1[nH]c2,hydroxytryptophan
#CC(=CCC/C(=C/CC/C(=C/CC/C=C(/CC/C=C(/CCC=C(C)C)\C)\C)/C)/C)C,squalene
```

Por razões didáticas, solventes e solutos estão separados em arquivos distintos.
Boas práticas de programação e design de protocolo, entretanto, dão preferência
a formas mais compactas de apresentação dos dados. O usuário, se quiser, pode
modificar os scripts e deixá-los mais ao gosto do cliente.

## Organização dos arquivos

O script `organize.sh` gerará três pastas:

- `001.solutes`
- `002.solvents`
- `003.initial_boxes`

As pastas `001.solutes` e `002.solvents` apenas contém as estruturas e parâmetros
de cada molécula individualmente. Os arquivos de input para as simulações de
dinâmica molecular estão em `003.initial_boxes` em diretórios nomeados como
`solute_in_solvent`:

```
biphenyl_in_acetaldehyde         
biphenyl_in_benzene              
biphenyl_in_chloroform           
biphenyl_in_ethanol              
biphenyl_in_methanol             
biphenyl_in_tetrahydrofuran      
```

E dentro de cada uma dessas pastas você encontrará entre os arquivos:
```
solute_in_solvent.gro
solute_in_solvent.top
minimization.mdp
equil_nvt.mdp
equil_npt.mdp
equil_npt2.mdp
prod.mdp
```
Esses são os arquivos importantes para as simulações: `solute_in_solvent.gro`
contém as coordenadas iniciais do sistema; `solute_in_solvent.top` contém os
parâmetros do campo de forças de cada molécula. Os arquivos `.mdp` contém os
parâmetros das simulações de dinâmica molecular e merecem uma explicação à parte.

## Rodando uma simulação de Dinâmica Molecular

São três os estágios para a produção de uma boa trajetória de dinâmica
molecular:

### Minimização
A energia em função das coordenadas atômicas no sistema é
minimizada de forma a reduzir efeitos energéticos indesejáveis de impedimentos
estéreos. Sem minimização, simulações de dinâmica molecular podem dar errado
("explodir"/"_blow up_") devido a forças elevadas que poderiam ter sido atenuadas
com a minimização.

Para fazer a minimização, digite no terminal:

```
gmx grompp -f minimization.mdp -c solute_in_solvent.gro -p solute_in_solvent.top -o minimization.tpr
gmx mdrun -deffnm minimization
```

O primeiro comando organiza a simulação, o segundo faz a dinâmica molecular ou
a minimização rodar no seu computador. A _flag_ `-f` indica o arquivo `.mdp`,
`-c` indica o arquivo `.gro` e `-p` indica o arquivo `.top`. O _output_ de
interesse para nós é `minimization.gro` (nomeado pela _flag_ `-deffnm`), que
contém a estrutura minimizada do sistema de interesse. `-o` é o _output_ com o
sistema pronto para o programa `mdrun`.

### Equilibração
Após a minimização a energia se encontra em um mínimo e o sistema
não necessariamente se encontra em uma configuração representativa. Além disso,
em dinâmica molecular, estamos interessados em propriedades de _ensemble_ e, mesmo
se a configuração inicial fosse quimicamente relevante, ela não teria o devido
valor estatístico.

São três os estágios de equilibração:
- __Equilibração de temperatura:__ `equil_nvt.mdp` define o estágio de equilibração
 da temperatura. Essa etapa é necessária porque a estrutura minimizada corresponderia
 ao sistema na temperatura de $0 \text{kelvin}$. Além disso, é no início desta etapa
 que velocidades são atribuídas a cada átomo, o que é imprescindível para a simulação.
 Para rodar esse estágio, lembre-se que a estrutura de partida está definida no
 _output_ da minimização:
 ```
 gmx grompp -f equil_nvt.mdp -c minimization.gro -p solute_in_solvent.top -o equil_nvt.tpr
 gmx mdrun -deffnm equil_nvt
 ```

- __Equilibração de pressão 1:__ `equil_npt.mdp` define o primeiro estágio de
equilibração da pressão. São necessários dois estágios porque o barostato de Berendsen
usado neste estágio, ajuda a caixa de simulação ficar com um volume adequado, mas
não produz um ensemble isotérmico-isobárico adequado. Para isso necessitamos de um
segundo estágio de equilibração da pressão. Para isso, usamos o arquivo `.gro`
produzido pela etapa de equilíbrio da temperatura, pois cada átomo além de conter
suas coordenadas atômicas, também tem as componentes de sua velocidade em uma
temperatura adequada:
```
gmx grompp -f equil_npt.mdp -c equil_nvt.gro -p solute_in_solvent.top -o equil_npt.tpr
gmx mdrun -deffnm equil_npt
```

- __Equilibração de pressão 2:__ `equil_npt2.mdp` define o segundo estágio de
equilibração da pressão. Neste estágio, usamos o barostato de Parrinello-Rahman,
que produz um ensemble isotérmico-isobárico adequado. A diferença entre o estágio
anterior e este é que o comprimento dos lados da caixa podem variar de forma
independente, enquanto isso não era o caso do estágio anterior. Para obter os
resultados:
```
gmx grompp -f equil_npt2.mdp -c equil_npt.gro -p solute_in_solvent.top -o equil_npt2.tpr
gmx mdrun -deffnm equil_npt2
```

### Produção

O estágio mais importante de uma simulação de dinâmica molecular é a amostragem
das configurações de equilíbrio do sistema. Isso é feito na etapa de produção.  
`prod.mdp` define este estágio. A produção é significantemente maior que os
demais passos; a amostragem de configurações de moléculas pequenas em água
costuma exigir 2500000 passos, implicando em um tempo computacional de algumas
horas.

A produção resulta na criação de uma trajetória contendo configurações no _ensemble_
isotérmico-isobárico (NPT) que podem ser usadas para calcular propriedades de
equilíbrio do sistema. A simulação é rodada com os seguintes comandos:
```
gmx grompp -f prod.mdp -c equil_npt2.gro -p solute_in_solvent.top -o prod.tpr
gmx mdrun -deffnm prod
```

Todos os dados da trajetória estão armazenados em arquivos `.trr` e
`.edr` que podem ser analisados com outros programas incluídos no GROMACS.

### Troubleshooting
Caso alguma etapa tenha falhado, devemos conferir o arquivo `.log` produzido.
