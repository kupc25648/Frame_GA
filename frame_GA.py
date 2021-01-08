'''
==================================================================
GA file for frame brance
==================================================================
# Prof Yamamoto and Prof Minakawa method
# 6bit represent [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
# But swap
# 3bit = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7]
# 1bit = [+0.0,+0.1]
# 1bit = [+0.0,+0.1]
# 1bit = [+0.0,+0.1]
'''
from FEM_frame import *
from frame_GEN import *
import numpy as np
import matplotlib.pyplot as plt

def plotforgame(data,pop,iter_num,init,final,detail=''):
    objective = detail+'Strain energy (Nm)'
    name_of_file = 'GA_population{}_iteration{}_Initial_{}_Final_{}_'.format(
        int(pop),int(iter_num),round(init,2),round(final,2))+detail+'.png'
    plt.ylabel('{}'.format(objective))
    plt.xlabel("step")
    plt.plot(data)
    plt.savefig(name_of_file)
    plt.close("all")


class Gene:
    def __init__(self,gen_model):
        self.name = 'GENE'
        self.code = []
        self.gen_model = gen_model
        self.model = gen_model.model
        self.model_gene = gen_model.model
        self.add_mark = 0

    def gen_gene(self):
        for i in range(len(self.model.nodes)):
            #self.code.append(np.random.choice([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]))


            #self.code.append(np.random.choice([0.0,0.2,0.4,0.6,0.8,1.0]))

            # genearate 6 bits
            bit = []
            for i in range(6):
                 bit.append(np.random.choice([0,1]))
            self.code.append(bit)


class Population:
    def __init__(self,objective,num_gene,nx,nz,span,diameter,lx,ly,lz,young,c1,c2,c3,c4,c5,c6,c7,forgame=None,game_max_y_val=None,brace=None):
        self.name = 'GENE'
        self.gene = []
        self.gene_check = []
        self.avg_u = 0
        self.zoo = []

        self.init = 0
        self.avg_init = 0
        self.avg_curr = 0
        self.objective = objective
        self.num_gene = num_gene

        # model parameter
        self.nx=nx
        self.nz=nz
        self.span=span
        self.diameter=diameter
        self.lx=lx
        self.ly=ly
        self.lz=lz
        self.young=young
        self.c1=c1
        self.c2=c2
        self.c3=c3
        self.c4=c4
        self.c5=c5
        self.c6=c6
        self.c7=c7
        self.brace_type = None
        self.forgame = forgame
        self.game_max_y_val = 1
        self.dummy_str1 = None
        self.dummy_str2 = None

        self._addgene(self.num_gene)
        # initialize dummy structure
        self.dummy_str1 = gen_model(self.nx,self.nz,self.span,self.diameter,self.lx,self.ly,self.lz,self.young,self.c1,self.c2,self.c3,self.c4,self.c5,self.c6,self.c7,self.forgame,self.game_max_y_val,self.brace_type)
        self.dummy_str2 = gen_model(self.nx,self.nz,self.span,self.diameter,self.lx,self.ly,self.lz,self.young,self.c1,self.c2,self.c3,self.c4,self.c5,self.c6,self.c7,self.forgame,self.game_max_y_val,self.brace_type)


    def _addgene(self,num_gene):
        allzoo = []
        for i in range(num_gene):
            allzoo.append('zoo{}'.format(i+1))
            allzoo[-1] = gen_model(self.nx,self.nz,self.span,self.diameter,self.lx,self.ly,self.lz,self.young,self.c1,self.c2,self.c3,self.c4,self.c5,self.c6,self.c7,self.forgame,self.game_max_y_val,self.brace_type)

        all_gene = []
        for i in range(num_gene):
            all_gene.append('g{}'.format(i+1))
            all_gene[-1] = Gene(allzoo[i])
            all_gene[-1].gen_gene()
            self.gene.append(all_gene[-1])

    def genstr(self):
        #interpret gene

        self.gene_check = []
        for i in range(len(self.gene)):

            self.gene[i].model_gene.restore()
            mark = 0

            for j in range(len(self.gene[i].code)):


                # transform code(6bit) to z-val
                # transform 3bits
                holder = 0
                if (self.gene[i].code[j][0],self.gene[i].code[j][2],self.gene[i].code[j][4]) == (0,0,0):
                    holder = 0
                elif (self.gene[i].code[j][0],self.gene[i].code[j][2],self.gene[i].code[j][4]) == (0,0,1):
                    holder = 0.1
                elif (self.gene[i].code[j][0],self.gene[i].code[j][2],self.gene[i].code[j][4]) == (0,1,0):
                    holder = 0.2
                elif (self.gene[i].code[j][0],self.gene[i].code[j][2],self.gene[i].code[j][4]) == (1,0,0):
                    holder = 0.3
                elif (self.gene[i].code[j][0],self.gene[i].code[j][2],self.gene[i].code[j][4]) == (0,1,1):
                    holder = 0.4
                elif (self.gene[i].code[j][0],self.gene[i].code[j][2],self.gene[i].code[j][4]) == (1,0,1):
                    holder = 0.5
                elif (self.gene[i].code[j][0],self.gene[i].code[j][2],self.gene[i].code[j][4]) == (1,1,0):
                    holder = 0.6
                else:
                    holder = 0.7
                holder += 0.1*self.gene[i].code[j][1]
                holder += 0.1*self.gene[i].code[j][3]
                holder += 0.1*self.gene[i].code[j][5]

                self.gene[i].model.nodes[j].coord[1] = holder
                if self.gene[i].model.nodes[j].res[0] == 1:
                    self.gene[i].model.nodes[j].coord[1] = 0.0
                mark += 1
            self.gene[i].model_gene.gen_all()

            self.gene[i].add_mark = mark

            self.gene_check.append([self.gene[i],self.gene[i].model_gene.U_full[0][0]])
        '''
        if self.init == 0:
            sum_e = 0
            for i in range(len(self.gene_check)):
                sum_e += self.gene_check[i][1]
            self.avg_init = sum_e/len(self.gene_check)
            self.init = 1


        self.gene_check = sorted(self.gene_check, key = lambda x: x[1])
        '''
        #generate structure from gene

    def fitness(self,count):
        # Sort strain energy
        if self.init == 0:
            sum_e = 0
            for i in range(len(self.gene_check)):
                sum_e += self.gene_check[i][1]
            self.avg_init = sum_e/len(self.gene_check)
            self.init = 1
        else:
            sum_e = 0
            for i in range(len(self.gene_check)):
                sum_e += self.gene_check[i][1]

        # SORT low TO high
        self.gene_check = sorted(self.gene_check, key = lambda x: x[1])

        iter_min_e = self.gene_check[0][1]
        iter_avg_e = sum_e/len(self.gene_check)

        print('AVG INIT ENERGY {} || AVG ENERGY NOW {} || MIN ENERGY NOW {}'.format(self.avg_init,sum_e/len(self.gene_check),self.gene_check[0][1]))
        return iter_min_e,iter_avg_e
        #self.gene_check = self.gene_check[:-int(len(self.gene_check)*0.3)]

        #print('CURRENT FIT ENERGY {}'.format(self.avg_curr))
        #print('FIT NUM {}'.format(len(self.gene_check)))


class GA:
    def __init__(self,population,pop_multate_rate,gene_multate_rate,cross_rate,max_pop,survive_rate,iter = 100):
        self.pool = population
        self.pool_backup = None
        self.con =0
        self.count = iter
        self.pop_multate_rate = pop_multate_rate
        self.gene_multate_rate = gene_multate_rate
        self.cross_rate = cross_rate
        self.max_pop = max_pop
        self.survive_rate = survive_rate

    def selection(self):
        self.pool.genstr()
        data1,data2 = self.pool.fitness(self.count)

        return data1,data2

    def cross(self):
        '''
        --------------------
        Algorithm
        --------------------
        1. Sort Best -> Worst (already done)
        2. Pair [1,2], [3,4], ...
        3. Crossover from pairs get [1',2'],[3',4'],...
        4. Mutation (each gene hase mutation rate of 1/len(chromosome) in this case 1/6)
        5. Group[1,2,1',2'],[3,4,3',4'],...
        6. Sort[1,2,1',2'],[3,4,3',4'],...
        7. Kill 2 worst in groups
        8. iteration over
        '''
        # Pair
        group = []
        for i in range(0,int(self.max_pop),2):
            # [code, energy]
            group.append([[self.pool.gene_check[i][0].code,self.pool.gene_check[i][1]],[self.pool.gene_check[i+1][0].code,self.pool.gene_check[i+1][1]]])
        # Crossover Mutation and make group
        for i in range(len(group)):
            # Crossover

            off1 = group[i][0][0].copy()
            off2 = group[i][1][0].copy()


            '''
            # type 2
            #--------------------------
            off1[12:] = group[i][1][0][12:]
            off2[12:] = group[i][0][0][12:]
            for j in range(len(group[i][0])):
            #--------------------------
            '''
            # type 1
            #--------------------------
            for j in range(len(group[i][0])):
                off1[j][3:] = group[i][1][0][j][3:]
                off2[j][3:] = group[i][0][0][j][3:]
            #--------------------------
                # Mutation
                for k in range(len(off1[j])):
                    rate1 = random.random()
                    rate2 = random.random()
                    if rate1 <= 1/6:
                        if off1[j][k] == 0:
                            off1[j][k] = 1
                        else:
                            off1[j][k] = 0
                    if rate2 <= 1/6:
                        if off2[j][k] == 0:
                            off2[j][k] = 1
                        else:
                            off2[j][k] = 0

            self.pool.dummy_str1.model.restore()
            self.pool.dummy_str2.model.restore()
            for j in range(len(off1[j])):
                # off1
                # transform code(6bit) to z-val
                # transform 3bits
                holder1 = 0
                if (off1[j][0],off1[j][2],off1[j][4]) == (0,0,0):
                    holder1 = 0
                elif (off1[j][0],off1[j][2],off1[j][4]) == (0,0,1):
                    holder1 = 0.1
                elif (off1[j][0],off1[j][2],off1[j][4]) == (0,1,0):
                    holder1 = 0.2
                elif (off1[j][0],off1[j][2],off1[j][4]) == (1,0,0):
                    holder1 = 0.3
                elif (off1[j][0],off1[j][2],off1[j][4]) == (0,1,1):
                    holder1 = 0.4
                elif (off1[j][0],off1[j][2],off1[j][4]) == (1,0,1):
                    holder1 = 0.5
                elif (off1[j][0],off1[j][2],off1[j][4]) == (1,1,0):
                    holder1 = 0.6
                else:
                    holder1 = 0.7
                holder1 += 0.1*off1[j][1]
                holder1 += 0.1*off1[j][3]
                holder1 += 0.1*off1[j][5]

                self.pool.dummy_str1.model.nodes[j].coord[1] = holder1
                if self.pool.dummy_str1.model.nodes[j].res[0] == 1:
                    self.pool.dummy_str1.model.nodes[j].coord[1] = 0.0

                # off2
                # transform code(6bit) to z-val
                # transform 3bits
                holder2 = 0
                if (off2[j][0],off2[j][2],off2[j][4]) == (0,0,0):
                    holder2 = 0
                elif (off2[j][0],off2[j][2],off2[j][4]) == (0,0,1):
                    holder = 0.1
                elif (off2[j][0],off2[j][2],off2[j][4]) == (0,1,0):
                    holder2 = 0.2
                elif (off2[j][0],off2[j][2],off2[j][4]) == (1,0,0):
                    holder2 = 0.3
                elif (off2[j][0],off2[j][2],off2[j][4]) == (0,1,1):
                    holder2 = 0.4
                elif (off2[j][0],off2[j][2],off2[j][4]) == (1,0,1):
                    holder2 = 0.5
                elif (off2[j][0],off2[j][2],off2[j][4]) == (1,1,0):
                    holder2 = 0.6
                else:
                    holder2 = 0.7
                holder2 += 0.1*off1[j][1]
                holder2 += 0.1*off1[j][3]
                holder2 += 0.1*off1[j][5]

                self.pool.dummy_str2.model.nodes[j].coord[1] = holder2
                if self.pool.dummy_str2.model.nodes[j].res[0] == 1:
                    self.pool.dummy_str2.model.nodes[j].coord[1] = 0.0

            self.pool.dummy_str1.model.gen_all()
            self.pool.dummy_str2.model.gen_all()

            # Group
            group[i].append([off1,self.pool.dummy_str1.model.U_full[0][0]])
            group[i].append([off2,self.pool.dummy_str2.model.U_full[0][0]])

            # Sort Group
            group[i] = sorted(group[i], key = lambda x: x[1])

        self.pool.gene = []

        for i in range(len(group)):
            # first place
            self.pool._addgene(1)
            self.pool.gene[-1].code = group[i][0][0]
            # second place
            self.pool._addgene(1)
            self.pool.gene[-1].code = group[i][1][0]

        out_num = len(self.pool.gene)
        return out_num

    # put before cross
    def cut(self,percent):
        #dummy_gene_check = self.pool.gene_check.copy()
        #dummy_gene_check.reverse()
        #cut = int((1-percent)*len(self.pool.gene_check))

        #self.pool.gene_check[cut:] = dummy_gene_check[cut:]
        pass


    def check_con(self):

        if (sum(self.pool.gene_check)/len(self.pool.gene_check)) <= ((1-self.pool.objective)*self.pool.avg_init):
            print('FOUND SOLUTION')
            for i in range(len(self.pool.gene)):
                if self.pool.gene[i].gen_model.model.U_full[0][0] <= (sum(self.pool.gene_check)/len(self.pool.gene_check)):
                    print(self.pool.gene[i].code)
                    self.pool.gene[i].gen_model.model =  self.pool.gene[i].model_gene
                    self.pool.gene[i].gen_model.render('Optimum solution',self.count,self.pool.gene[i].gen_model.model.U_full[0][0])




    def run(self):
        #self.count = 1
        num  = 0
        min_story = []
        avg_story = []
        while self.count != num:
            print('ITERATION {} NUM_GENE {}'.format(num+1,len(self.pool.gene)))
            data1,data2 = self.selection()
            min_story.append(data1)
            avg_story.append(data2)
            if num%5 == 3:
                self.cut(0.3)
            out = self.cross()
            #self.check_con()
            num+=1
            print('---------------------------------------------------')

        plotforgame(min_story,self.max_pop,self.count,min_story[0],min_story[-1],'MIN')
        plotforgame(avg_story,self.max_pop,self.count,min_story[0],min_story[-1],'AVG')

        for i in range(len(self.pool.gene)):

            # some criteria

            for i in range(len(self.pool.gene)):
                self.pool.gene[i].model_gene.restore()
                for j in range(len(self.pool.gene[i].code)):
                    holder = 0
                    if (self.pool.gene[i].code[j][0],self.pool.gene[i].code[j][2],self.pool.gene[i].code[j][4]) == (0,0,0):
                        holder = 0
                    elif (self.pool.gene[i].code[j][0],self.pool.gene[i].code[j][2],self.pool.gene[i].code[j][4]) == (0,0,1):
                        holder = 0.1
                    elif (self.pool.gene[i].code[j][0],self.pool.gene[i].code[j][2],self.pool.gene[i].code[j][4]) == (0,1,0):
                        holder = 0.2
                    elif (self.pool.gene[i].code[j][0],self.pool.gene[i].code[j][2],self.pool.gene[i].code[j][4]) == (1,0,0):
                        holder = 0.3
                    elif (self.pool.gene[i].code[j][0],self.pool.gene[i].code[j][2],self.pool.gene[i].code[j][4]) == (0,1,1):
                        holder = 0.4
                    elif (self.pool.gene[i].code[j][0],self.pool.gene[i].code[j][2],self.pool.gene[i].code[j][4]) == (1,0,1):
                        holder = 0.5
                    elif (self.pool.gene[i].code[j][0],self.pool.gene[i].code[j][2],self.pool.gene[i].code[j][4]) == (1,1,0):
                        holder = 0.6
                    else:
                        holder = 0.7
                    holder += 0.1*self.pool.gene[i].code[j][1]
                    holder += 0.1*self.pool.gene[i].code[j][3]
                    holder += 0.1*self.pool.gene[i].code[j][5]

                    self.pool.gene[i].model.nodes[j].coord[1] = holder



                    #self.pool.gene[i].model.nodes[j].coord[1] = self.pool.gene[i].code[j]



                    if self.pool.gene[i].model.nodes[j].res[0] == 1:
                        self.pool.gene[i].model.nodes[j].coord[1] = 0.0
                self.pool.gene[i].model_gene.gen_all()
            min_e = []
            for i in range(len(self.pool.gene)):
                min_e.append(self.pool.gene[i].gen_model.model.U_full[0][0])
            min_e = min(min_e)

            for i in range(len(self.pool.gene)):
                if self.pool.gene[i].gen_model.model.U_full[0][0] <= min_e*1.05:
                    print(self.pool.gene[i].code)
                    self.pool.gene[i].gen_model.render('Optimum solution',self.count,self.pool.gene[i].gen_model.model.U_full[0][0])






