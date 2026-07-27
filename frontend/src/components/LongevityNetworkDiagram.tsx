import { useMemo } from 'react';
import ReactFlow, {
  Background,
  BackgroundVariant,
  Controls,
  Handle,
  MiniMap,
  Position,
  ReactFlowProvider,
  useEdgesState,
  useNodesState,
  type Edge,
  type Node,
  type NodeProps,
} from 'reactflow';
import 'reactflow/dist/style.css';
import {
  Activity,
  Brain,
  Dna,
  FlaskConical,
  GitBranch,
  Sparkles,
} from 'lucide-react';

export type LongevityNodeCategory =
  | 'core'
  | 'pathophysiology'
  | 'intersection'
  | 'mechanism'
  | 'axis'
  | 'genes'
  | 'outcome';

export type LongevityNodeData = {
  title: string;
  subtitle?: string;
  category: LongevityNodeCategory;
};

const categoryClass: Record<LongevityNodeCategory, string> = {
  core: 'border-l-blue-500 bg-blue-50',
  pathophysiology: 'border-l-blue-500 bg-blue-50',
  intersection: 'border-l-blue-500 bg-blue-50',
  mechanism: 'border-l-amber-500 bg-amber-50',
  axis: 'border-l-emerald-500 bg-emerald-50',
  genes: 'border-l-purple-500 bg-purple-50',
  outcome: 'border-l-emerald-500 bg-emerald-50',
};

const categoryIcon: Record<LongevityNodeCategory, typeof Sparkles> = {
  core: Sparkles,
  pathophysiology: Brain,
  intersection: GitBranch,
  mechanism: FlaskConical,
  axis: Activity,
  genes: Dna,
  outcome: Brain,
};

function CustomNode({ data, selected }: NodeProps<LongevityNodeData>) {
  const Icon = categoryIcon[data.category];
  const border = categoryClass[data.category];

  return (
    <div
      className={[
        'relative min-w-[200px] max-w-[280px] rounded-lg border border-slate-200/80 bg-white pl-3 pr-3 pb-3 pt-3',
        'shadow-sm transition-all duration-200 ease-out',
        'hover:z-10 hover:scale-[1.02] hover:shadow-md',
        selected ? 'ring-2 ring-slate-400/40 ring-offset-2' : '',
        border,
      ].join(' ')}
    >
      <Handle
        type="target"
        position={Position.Top}
        className="h-2 w-2 border border-slate-300 bg-white"
      />
      <div className="flex items-start justify-between gap-2">
        <div className="min-w-0 flex-1">
          <p className="text-sm font-semibold leading-snug text-slate-800">
            {data.title}
          </p>
          {data.subtitle ? (
            <p className="mt-1 text-xs leading-relaxed text-slate-600">
              {data.subtitle}
            </p>
          ) : null}
        </div>
        <Icon
          className="mt-0.5 h-4 w-4 shrink-0 text-slate-400"
          strokeWidth={1.75}
          aria-hidden
        />
      </div>
      <Handle
        type="source"
        position={Position.Bottom}
        className="h-2 w-2 border border-slate-300 bg-white"
      />
    </div>
  );
}

const nodeTypes = { longevityCustom: CustomNode };

function buildGraph(): { nodes: Node<LongevityNodeData>[]; edges: Edge[] } {
  const nodes: Node<LongevityNodeData>[] = [
    {
      id: 'core',
      type: 'longevityCustom',
      position: { x: 380, y: 0 },
      data: {
        category: 'core',
        title: 'Exceptional Longevity: Active Genetic Resilience',
      },
    },
    {
      id: 'adpd',
      type: 'longevityCustom',
      position: { x: 340, y: 110 },
      data: {
        category: 'pathophysiology',
        title: 'AD & PD Pathophysiology',
      },
    },
    {
      id: 'intersection',
      type: 'longevityCustom',
      position: { x: 260, y: 220 },
      data: {
        category: 'intersection',
        title: 'Dysregulated Network: 41 Genes, 70 LA-SNPs',
      },
    },
    {
      id: 'qualitative',
      type: 'longevityCustom',
      position: { x: 120, y: 360 },
      data: {
        category: 'mechanism',
        title: 'Qualitative Effects (Protein Efficiency & Stability)',
        subtitle: 'Primary mechanistic path',
      },
    },
    {
      id: 'eqtl',
      type: 'longevityCustom',
      position: { x: 560, y: 360 },
      data: {
        category: 'mechanism',
        title: 'eQTL Effects (Gene Expression)',
        subtitle: 'Secondary mechanistic path',
      },
    },
    {
      id: 'axis1',
      type: 'longevityCustom',
      position: { x: 0, y: 520 },
      data: { category: 'axis', title: 'Protein Homeostasis' },
    },
    {
      id: 'genes1',
      type: 'longevityCustom',
      position: { x: 0, y: 640 },
      data: {
        category: 'genes',
        title: 'Target genes',
        subtitle: 'HSPA1A, HSPA1B, HSPA1L',
      },
    },
    {
      id: 'outcome1',
      type: 'longevityCustom',
      position: { x: 0, y: 760 },
      data: {
        category: 'outcome',
        title: 'Outcome',
        subtitle: 'Proteotoxic stress resilience',
      },
    },
    {
      id: 'axis2',
      type: 'longevityCustom',
      position: { x: 280, y: 520 },
      data: { category: 'axis', title: 'Lipid Metabolism & Mitochondria' },
    },
    {
      id: 'genes2',
      type: 'longevityCustom',
      position: { x: 280, y: 640 },
      data: {
        category: 'genes',
        title: 'Target genes',
        subtitle: 'CETP (rs5882), NDUFS1',
      },
    },
    {
      id: 'outcome2',
      type: 'longevityCustom',
      position: { x: 280, y: 760 },
      data: {
        category: 'outcome',
        title: 'Outcome',
        subtitle: 'Altered HDL & limits ROS',
      },
    },
    {
      id: 'axis3',
      type: 'longevityCustom',
      position: { x: 560, y: 520 },
      data: { category: 'axis', title: 'Immune Regulation' },
    },
    {
      id: 'genes3',
      type: 'longevityCustom',
      position: { x: 560, y: 640 },
      data: {
        category: 'genes',
        title: 'Target genes',
        subtitle: 'NLRC5, HLA-DQB1',
      },
    },
    {
      id: 'outcome3',
      type: 'longevityCustom',
      position: { x: 560, y: 760 },
      data: {
        category: 'outcome',
        title: 'Outcome',
        subtitle: "Attenuates 'Inflammaging'",
      },
    },
    {
      id: 'axis4',
      type: 'longevityCustom',
      position: { x: 840, y: 520 },
      data: { category: 'axis', title: 'Antioxidant Factors' },
    },
    {
      id: 'genes4',
      type: 'longevityCustom',
      position: { x: 840, y: 640 },
      data: {
        category: 'genes',
        title: 'Target genes',
        subtitle: 'HMOX1, GPX1',
      },
    },
    {
      id: 'outcome4',
      type: 'longevityCustom',
      position: { x: 840, y: 760 },
      data: {
        category: 'outcome',
        title: 'Outcome',
        subtitle: 'Preserves mitochondrial structure',
      },
    },
  ];

  const edges: Edge[] = [
    { id: 'e-core-adpd', source: 'core', target: 'adpd', animated: false },
    { id: 'e-adpd-intersection', source: 'adpd', target: 'intersection' },
    {
      id: 'e-intersection-qualitative',
      source: 'intersection',
      target: 'qualitative',
      style: { strokeWidth: 2 },
    },
    {
      id: 'e-intersection-eqtl',
      source: 'intersection',
      target: 'eqtl',
      style: { strokeWidth: 1.5, strokeDasharray: '6 4' },
    },
    { id: 'e-qual-axis1', source: 'qualitative', target: 'axis1' },
    { id: 'e-qual-axis2', source: 'qualitative', target: 'axis2' },
    { id: 'e-qual-axis3', source: 'qualitative', target: 'axis3' },
    { id: 'e-qual-axis4', source: 'qualitative', target: 'axis4' },
    { id: 'e-a1-g1', source: 'axis1', target: 'genes1' },
    { id: 'e-g1-o1', source: 'genes1', target: 'outcome1' },
    { id: 'e-a2-g2', source: 'axis2', target: 'genes2' },
    { id: 'e-g2-o2', source: 'genes2', target: 'outcome2' },
    { id: 'e-a3-g3', source: 'axis3', target: 'genes3' },
    { id: 'e-g3-o3', source: 'genes3', target: 'outcome3' },
    { id: 'e-a4-g4', source: 'axis4', target: 'genes4' },
    { id: 'e-g4-o4', source: 'genes4', target: 'outcome4' },
  ];

  return { nodes, edges };
}

function LongevityNetworkDiagramInner() {
  const { nodes: initialNodes, edges: initialEdges } = useMemo(
    () => buildGraph(),
    [],
  );
  const [nodes, , onNodesChange] = useNodesState(initialNodes);
  const [edges, , onEdgesChange] = useEdgesState(initialEdges);

  const minimapNodeColor = (node: Node<LongevityNodeData>) => {
    switch (node.data.category) {
      case 'core':
      case 'pathophysiology':
      case 'intersection':
        return '#3b82f6';
      case 'mechanism':
        return '#f59e0b';
      case 'genes':
        return '#a855f7';
      default:
        return '#10b981';
    }
  };

  return (
    <div
      data-longevity-diagram="capture"
      className="relative h-full min-h-[720px] w-full rounded-xl border border-slate-200 bg-slate-50/80"
    >
      <ReactFlow
        nodes={nodes}
        edges={edges}
        onNodesChange={onNodesChange}
        onEdgesChange={onEdgesChange}
        nodeTypes={nodeTypes}
        fitView
        fitViewOptions={{ padding: 0.15 }}
        minZoom={0.35}
        maxZoom={1.5}
        className="rounded-xl"
        defaultEdgeOptions={{
          style: { stroke: '#64748b', strokeWidth: 1.25 },
        }}
      >
        <Background
          variant={BackgroundVariant.Dots}
          gap={18}
          size={1}
          color="#cbd5e1"
        />
        <Controls
          className="overflow-hidden rounded-md border border-slate-200 bg-white shadow-sm"
          showInteractive={false}
        />
        <MiniMap
          nodeColor={minimapNodeColor}
          maskColor="rgb(248 250 252 / 0.85)"
          className="rounded-md border border-slate-200 bg-white"
        />
      </ReactFlow>
      <div className="pointer-events-none absolute bottom-3 left-3 flex flex-wrap gap-3 rounded-md border border-slate-200 bg-white/95 px-3 py-2 text-[10px] text-slate-600 shadow-sm backdrop-blur-sm">
        <span className="flex items-center gap-1.5">
          <span className="h-2 w-2 rounded-sm bg-blue-500" />
          Core / intersection
        </span>
        <span className="flex items-center gap-1.5">
          <span className="h-2 w-2 rounded-sm bg-amber-500" />
          Mechanism
        </span>
        <span className="flex items-center gap-1.5">
          <span className="h-2 w-2 rounded-sm bg-emerald-500" />
          Axes / outcomes
        </span>
        <span className="flex items-center gap-1.5">
          <span className="h-2 w-2 rounded-sm bg-purple-500" />
          Target genes
        </span>
        <span className="flex items-center gap-1.5 border-l border-slate-200 pl-3">
          <span className="font-medium text-slate-700">Solid</span>
          qualitative path
        </span>
        <span className="flex items-center gap-1.5">
          <svg width={28} height={6} aria-hidden className="text-slate-500">
            <line
              x1={0}
              y1={3}
              x2={28}
              y2={3}
              stroke="currentColor"
              strokeWidth={1.5}
              strokeDasharray="4 3"
            />
          </svg>
          eQTL path
        </span>
      </div>
    </div>
  );
}

export function LongevityNetworkDiagram() {
  return (
    <ReactFlowProvider>
      <div className="relative h-full w-full">
        <LongevityNetworkDiagramInner />
      </div>
    </ReactFlowProvider>
  );
}

export default LongevityNetworkDiagram;
