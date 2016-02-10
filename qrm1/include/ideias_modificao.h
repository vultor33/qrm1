/*

WARNING!!!!!
A FUNÇÃO Q-MONOPOLO LEVA A UM ZERO DIFERENTE DO 
.


OVERVIEW DO PROGRAMA

- Objeto global que seta os parâmetros* (PARAM)

- Existe um objeto que lê a molecula (Molecule)

- ScfCycle e o objeto que faz o ciclo SCF
  se você coloca o numero 1 ele faz a inicializacao
  calcula as integrais, chama o objeto Core_Matrix_calculation
  e etc.

- ScfProcedure controla tudo.

* WARNING!!! Repulsao de nucleos, caso seja menor que zero, retorna zero.

* WARNING!!! PARAM com gerador de ro1 e ro2 errados.

* WARNING!!! Pseudo diagonalizacao com valores ligeiramente
             diferentes.

Planejamento

- Criar um input que escolha o atomo a ser parametrizado. (no main)
- Automatizar a transferência de qRM1 -> para RM1 se o q inserido
for um.
  [Falta só testar]

- UHF
  (Criar um objeto abstrato que faz os métodos UHF, RHF e etc. se tiver outros).




- Preciso sinalizar que a leitura do input foi concluida com sucesso.







- Testar repulsao de core no molds.
   - N consigo entender essas coordenadas internas do molds
     desisto temporariariamente.

- Arrumar a otimizazao
 
- Transformar PARAM em um objeto normal.

- Montar um esquema para verificar as equações de meu programa.

IDEIAS

1. Nas classe Diagonaliza_Eigen e Diagonaliza_Matriz
   poderia ser adicionado uns inline getautovalor()
   para reduzir os membros publicos
  
2. Nos modulos de diagonalizao podia ter bem mais
	referências.

4. Estamos considerando apenas 4 orbitais para todos
   os atomos, os tipos tambem sao fixos,
   orbital 0 - S.
   orbital 1 - Px.
   orbital 2 - Py.
   orbital 3 - Pz.

6. Todos os parametros ja foram cadastrados em ua.
   for god sake!
   o alfa foi deixado em A-1 por enquanto.

8. Em atomos estao os parametros
   Em integrais_de_4_centros estao essas integrais
   E em Matriz_de_Core estao as integrais de 2 centros

9. ATENCAO -- ATENCAO!!!
   Nos operadores de troca e coulomb eu nao sei se e:
   lambda> sigma ou se vai pra todo lado. Coloquei
   pra todo lado.
   Isso acontece muito. Tipo, Gmilambda soma todos
   os ni sigma, nao sei se inclue os proprios milambda.

11. Lembra que o Mopac tem um shift na matriz de fock
    para acelerar a convergencia do scf.

14. ATENÇÃO!!! A conversao de coordenadas de internas para cartesianas
    ainda está com problemas.

16. Mais de uma molecula vai precisar de habilitar as duas integais
    de overlap. A para tau = 0 e a outra.
	
18. ATENCAO!!! As integrais de 4 centros estao com as mensagens de erro
    fechadas (estao return 0.0). Depois tem que trocar para exit.
	
19. O Stewart fala que com 5 iteracoes ele consegue atingir a convergencia
    do ro_2. Eu precisei de mais de 30.

20. Existe uma hora que a integral de 4 centros é calculada de forma invertida
    de acordo com o paper tem que trocar o sinal, mas na matriz do stewart
	nao acontece troca nenhuma. Avaliar isso com calma depois.

21. MAIS UM MENOS COM PROBLEMAS. De acordo com o artigo de integral de de 2 centros
    e necessario adicionar um menos no colchete quando ele é invertido. No MOPAC
	parece que essa regra não é seguida.
	<<<<STEWART>>>>
	Vou assumir que o artigo está errado.
	--> Eu acho que quando eu troco d1a - d1b eu ja troco o sinal tambem.

22. ATENCAOOOO!!! Troquei a ressonancia na mao mesmo. Cuidado
    --> parece que nao comuta. Um  lado e + e o outro e menos.
	
23. Tem umas razoes nas integrais de ressonancia 1/9 que e muito
    melhor colocar o numero ne. constexpr se quiser

24. o valor da ressonancia muda. s2px da 23 e o contrario -23
    verificar isso.

25. CARAI!!! Desconfiar dos sinais do spxpipi foi demais.

26. Em four center integrals:  
    // WARNING!!! Cases 14,15,16,17 had to be plus instead of minus.

27. Faz todo sentido que a atração nuclear tenha uma direção preferencial.
    Faz sentido que o cara do meio seja zero.
	Agora não achei nas equações:
	+ZB(spx|ss) para A>B
	-ZB(spx|ss) para A<B

28. Tem um ofstream no ciclo scf que poderia ser retirado.

29. AEEEEEEE!!! Lógico que muda o sinal do overlap!
    Tem o e(-imfi) nos harmonicos esfericos.
	Não tem que comutar mesmo naum, gg toboco.

30. Em vector<vector<double>> Diatomic_Overlaps::calc_overlap_AB(const Molecule &mol, int A, int B)
    tem um size=4 fixo é o numero de bases da rotacao.
	Cuidado.
	Ele aparece em varios outros lugares

31. IMPORTANTE!!!!
    Ele esta varrendo todas as integrais de 4 centros
	mesmo para A=B. Aí coloquei um return 0 em
	four center integrals na bagunca. Cuidado.

32. 2pz-2pz saida do kapa quando tau =0

33. Em Diatomic_Overlps --> esta A!=B
    mas acho que nao tem problema se for A<B

34. Em Matrix rotation esta definido duas matrizes
    rotacao, o trem mais feio que ja vi na minha vida.
	elas aparecem em diatomic_overlaps e em fourcenter
	(algum lugar).

35. Minha operacao de rotacao do overlap tem uma transposta
    la nada haver, diferente do pople. Mas é o que deu certo.
	pqp.

36. ATENCAO!!! A rotina que le o mop coloca a quantidade de
    atomos como a quantidade de linhas depois das 3 primeiras.
	entao não pode haver linha em branco quando a leitura acabar.

37. PROBLEMA CERTO!!!
    No RHF a integral mi lambda e de troca, no UHF esta estranho
	parece que mudou a integral, aff.

38. PROBLEMA CERTO!!!
    No UHF mi ni eu adicionei na mao uma soma dos outros atomos,
	o mesmo termo do mimi, como no RHF

39. DIIS como quase tudo esta com copia de matrizes desnecessarias.

40. DlibLinearSystem vai resolvero problema diretamente, seria
	bom criar um objeto de interface.

41. LENTIDAO: O determinante da matriz M do diis
              e calculada antes de fazer a inversa
			  para verificar se a matriz e singular.

			  
42. RETIRADO:  rhfMatrix_.electronic_energy(coreFockMatrix);
               de ScfCycle			  
			  
*/








