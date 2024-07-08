import pandas as pd
import math
import matplotlib.pyplot as plt




#df1 = pd.read_csv('/home/erik98m/masteroppgave/10agents_executions.csv', on_bad_lines='skip')
#df2 = pd.read_csv('/home/erik98m/alloc/n10.csv', on_bad_lines='skip')


def find_score(file1, file2):

    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    num_rows = min(len(df1), len(df2))
    l = []
    x_HC = []
    x_SA = []
    x_TS = []
    x_MS = []

    for i in range(num_rows):
        try:
            val1 = df1.at[i, 'HC']
            val2 = df1.at[i, 'SA']
            val3 = df1.at[i, 'TS']
            val5 = df1.at[i, 'TS milli']

            val4 = df2.at[i,'Result']

            hc = math.floor(val1) / val4
            sa = math.floor(val2) / val4
            ts = math.floor(val3) / val4
            ms = math.floor(val5) / val4

            x_HC.append(hc)
            x_SA.append(sa)
            x_TS.append(ts)
            x_MS.append(ms)

        except Exception as e:
            print(f"Error at index {i}: {e}")

    l.append(sum(x_HC)/100)
    l.append(sum(x_SA)/100)
    l.append(sum(x_TS)/100)
    l.append(sum(x_MS)/100)

    return l

def ga_score(file1, file2):

    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    num_rows = min(len(df1), len(df2))
    l= []
    x_GA = []

    for i in range(num_rows):
        try:
            val1 = df1.at[i, 'GA Solution']
            val4 = df2.at[i,'Result']

            ga = math.floor(val1) / val4
            x_GA.append(ga)

        except Exception as e:
            print(f"Error at index {i}: {e}")

    l.append(sum(x_GA)/100)

    return l

def ga_scoren2 (file1, file2):

    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    num_rows = min(len(df1), len(df2))
    l = []
    x_GA = []

    for i in range(num_rows):
        try:
            if df2.at[i,'Result'] == 0.0:
                x_GA.append(0.0)
               
            else:
                val1 = df1.at[i, 'GA Solution']
                val4 = df2.at[i,'Result']

                ga = math.floor(val1) / val4
                x_GA.append(ga)


        except Exception as e:
            print(f"Error at index {i}: {e}")

    l.append(sum(x_GA)/100)
  

    return l

def find_scoren2 (file1, file2):

    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    num_rows = min(len(df1), len(df2))
    l = []
    x_HC = []
    x_SA = []
    x_TS = []
    #x_MS = []

    for i in range(num_rows):
        try:
            if df2.at[i,'Result'] == 0.0:
                x_HC.append(0.0)
                x_SA.append(0.0)
                x_TS.append(0.0)
                #x_MS.append(0.0)
            else:
                val1 = df1.at[i, 'HC']
                val2 = df1.at[i, 'HC01']
                val3 = df1.at[i, 'SA100']
                #val5 = df1.at[i, 'SA5']

                val4 = df2.at[i,'Result']

                hc = math.floor(val1) / val4
                sa = math.floor(val2) / val4
                ts = math.floor(val3) / val4
                #ms = math.floor(val5) / val4

                x_HC.append(hc)
                x_SA.append(sa)
                x_TS.append(ts)
                #x_MS.append(ms)

        except Exception as e:
            print(f"Error at index {i}: {e}")

    l.append(sum(x_HC)/100)
    l.append(sum(x_SA)/100)
    l.append(sum(x_TS)/100)
    #l.append(sum(x_MS)/100)

    return l

# Compute the scores for each number of agents

results = [
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/2agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n2.csv'),
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/3agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n3.csv'),
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/4agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n4.csv'),
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/5agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n5.csv'),
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/6agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n6.csv'),
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/7agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n7.csv'),
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/8agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n8.csv'),
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/9agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n9.csv'),
    find_scoren2('/home/erik98m/masteroppgave/TS/eqmaxi/executions/10agents_executions.csv', '/home/erik98m/alloc/pythonPlots/TS/eqmaxi/n10.csv')
]

#results = [
    #ga_scoren2('/home/erik98m/masteroppgave/mmgacsv/2agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n2.csv'),
    #ga_score('/home/erik98m/masteroppgave/mmgacsv/3agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n3.csv'),
    #ga_score('/home/erik98m/masteroppgave/mmgacsv/4agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n4.csv'),
    #ga_score('/home/erik98m/masteroppgave/mmgacsv/5agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n5.csv'),
    #ga_score('/home/erik98m/masteroppgave/mmgacsv/6agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n6.csv'),
    #ga_score('/home/erik98m/masteroppgave/mmgacsv/7agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n7.csv'),
    #ga_score('/home/erik98m/masteroppgave/mmgacsv/8agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n8.csv'),
    #ga_score('/home/erik98m/masteroppgave/mmgacsv/9agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n9.csv'),
    #ga_score('/home/erik98m/masteroppgave/mmgacsv/10agents_executions.csv', '/home/erik98m/alloc/mmcsvga/n10.csv')
#]

print(results)

# Extract the number of agents
num_agents = list(range(2, 2 + len(results)))

# Extract individual lists for HC, SA, and TS solutions
#ga_scores = [result[0] for result in results]
hc_scores = [result[0] for result in results]
sa_scores = [result[1] for result in results]
ts_scores = [result[2] for result in results]
#ts_ms_scores = [result[3] for result in results]

# Plotting
plt.figure(figsize=(10, 6))
#plt.plot(num_agents, ga_scores, marker='o', label='GA Solution with maximin fitness')
plt.plot(num_agents, hc_scores, marker='o', label='TS 1s')
plt.plot(num_agents, sa_scores, marker='s', label='TS 1ms')
plt.plot(num_agents, ts_scores, marker='^', label='TS 0.1ms')
#plt.plot(num_agents, ts_ms_scores, marker='x', label='SA S5')

# Adding titles and labels
plt.title('Performance vs. Number of Agents')
plt.xlabel('Number of Agents')
plt.ylabel('Values')
plt.xticks(num_agents)
plt.ylim(0.0, max(max(hc_scores), max(sa_scores), max(ts_scores)) * 1.1)
#plt.ylim(0, max(ga_scores) * 1.1)
plt.legend()
plt.grid(True)
plt.show()    
