"""Input contracts that protect the nutrient-only comparison before execution."""
import copy
import json
import math
import random
import unittest

from .protocol import load_protocol, validate_protocol, digest, exact_steps
from .inputs import (make_pair_inputs, arm_inputs, assert_matched_arms, tick_schedule,
                     seed_streams, build_simulation, restore_founders)


class InputContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.protocol=load_protocol()
        cls.pair=make_pair_inputs(cls.protocol,11)

    def test_original_starvation_buffer_and_feasible_capacity_are_independent(self):
        p=self.protocol
        self.assertEqual(p['founders']['initial_food_pg'],2.5)
        self.assertTrue(all(s['food_capacity']==100 for s in p['species']))
        for f in self.pair['founders']:
            self.assertEqual(f['state']['food_level'],2.5)
            self.assertEqual(f['state']['food_capacity'],100)
        invalid=copy.deepcopy(p)
        invalid['species'][0]['food_capacity']=10
        with self.assertRaisesRegex(ValueError,'infeasible'):
            validate_protocol(invalid)

    def test_full_disc_geometry_and_random_labels_are_saved(self):
        pair=self.pair
        self.assertEqual(pair['geometry']['radius'],600)
        self.assertEqual(len(pair['founders']),500)
        self.assertEqual(sum(f['population']=='focal' for f in pair['founders']),5)
        focal=[f['agent_id'] for f in pair['founders'] if f['population']=='focal']
        self.assertNotEqual(focal,list(range(5)))
        for f in pair['founders']:
            state=f['state']
            self.assertLessEqual(math.hypot(state['x'],state['y']),600)
            self.assertEqual(state['base_size'],1)
            self.assertEqual(state['age'],0)
            self.assertIn('heading',state)
            self.assertIn('lifespan_variance',state)
            self.assertIn('state_timer',state)
            self.assertIn('q',state)
        # Area-uniform sample should not accidentally be uniform in radius.
        mean_r2=sum((f['state']['x']**2+f['state']['y']**2)/600**2 for f in pair['founders'])/500
        self.assertGreater(mean_r2,.4)
        self.assertLess(mean_r2,.6)

    def test_founder_generation_has_no_ecological_rng_side_effect(self):
        state=random.getstate()
        make_pair_inputs(self.protocol,4)
        self.assertEqual(random.getstate(),state)

    def test_species_container_order_does_not_construct_different_founders(self):
        reverse=copy.deepcopy(self.protocol)
        reverse['species'].reverse()
        self.assertEqual(make_pair_inputs(reverse,11),self.pair)

    def test_founders_are_independent_of_source_contrast(self):
        neutral=copy.deepcopy(self.protocol)
        neutral['feeding']['contrast_nodes']=0
        other=make_pair_inputs(neutral,11)
        self.assertEqual(other['founders'],self.pair['founders'])
        self.assertEqual(other['streams'],self.pair['streams'])

    def test_streams_are_named_distinct_saved_and_repeatable(self):
        a=seed_streams(15)
        self.assertEqual(a,seed_streams(15))
        self.assertEqual(len(set(v['seed'] for v in a.values())),len(a))
        for stream in a.values():
            self.assertIn('numpy_state',stream)
            self.assertIn('python_state',stream)
            self.assertIn('spawn_key',stream['seed_sequence'])
        self.assertNotEqual(a,seed_streams(16))
        with self.assertRaises(ValueError):seed_streams(True)

    def test_exact_feed_mass_times_and_physical_invariance(self):
        schedule=self.pair['feed_schedule']
        self.assertEqual(len(schedule),600)
        self.assertEqual([e['time_s'] for e in schedule],list(range(600)))
        self.assertEqual(sum(n['mass_pg'] for e in schedule for n in e['nodes']),15300)
        coarse=tick_schedule(schedule,.1);fine=tick_schedule(schedule,.05)
        self.assertEqual(list(coarse),[10*i for i in range(600)])
        self.assertEqual(list(fine),[20*i for i in range(600)])
        for i in range(600):self.assertEqual(coarse[10*i],fine[20*i])
        with self.assertRaisesRegex(ValueError,'exact multiple'):
            tick_schedule(schedule,.3)

    def test_each_batch_is_two_two_or_three_one_in_exact_same_circle(self):
        u,h=[arm_inputs(self.pair,x) for x in ('U','H')]
        self.assertTrue(assert_matched_arms(u,h))
        angle=u['orientation_rad'];nx,ny=math.cos(angle),math.sin(angle)
        for ue,he in zip(u['feed_schedule'],h['feed_schedule']):
            self.assertEqual(sum(n['x']*nx+n['y']*ny>=0 for n in ue['nodes']),2)
            self.assertEqual(sum(n['x']*nx+n['y']*ny>=0 for n in he['nodes']),3)
            self.assertEqual([n['node_id'] for n in ue['nodes']],[n['node_id'] for n in he['nodes']])
            for un,hn in zip(ue['nodes'],he['nodes']):
                self.assertEqual(un['mass_pg'],hn['mass_pg'])
                self.assertEqual(un['substrate_id'],hn['substrate_id'])
                self.assertAlmostEqual(math.hypot(un['x'],un['y']),math.hypot(hn['x'],hn['y']),places=10)

    def test_zero_contrast_has_identical_external_inputs_and_ecology(self):
        p=copy.deepcopy(self.protocol);p['feeding']['contrast_nodes']=0
        pair=make_pair_inputs(p,11);u,h=[arm_inputs(pair,x) for x in ('U','H')]
        self.assertTrue(assert_matched_arms(u,h))
        for k in ('arm','input_hash'):u.pop(k);h.pop(k)
        self.assertEqual(u,h)

    def test_between_arm_tampering_rejected_even_with_recomputed_hash(self):
        u,h=[arm_inputs(self.pair,x) for x in ('U','H')]
        mutations=[lambda a:a['feed_schedule'][1]['nodes'][0].update(mass_pg=10),
                   lambda a:a['founders'][0]['state'].update(food_level=80),
                   lambda a:a['protocol']['species'][0].update(speed=80),
                   lambda a:a['feed_schedule'][0]['nodes'][0].update(x=10),
                   lambda a:a['streams']['ecology'].update(seed=7)]
        for mutate in mutations:
            changed=copy.deepcopy(h);mutate(changed)
            changed['input_hash']=digest({k:v for k,v in changed.items() if k!='input_hash'})
            with self.subTest(mutate=mutate),self.assertRaises(ValueError):
                assert_matched_arms(u,changed)

    def test_unknown_and_clamped_parameters_are_rejected(self):
        cases=[('unknown',lambda p:p.update(magic=1)),
               ('species',lambda p:p['species'][0].update(greediness=1.2)),
               ('nan',lambda p:p['founders'].update(initial_food_pg=float('nan'))),
               ('clock',lambda p:p['clock'].update(dt_model_seconds=.03)),
               ('confound',lambda p:p['neutralized_radial_parameters'].update(MUCUS_MOTILITY_FACTOR=.5)),
               ('disabled',lambda p:p['species'][0].update(acetate_uptake=.2)),
               ('overlap',lambda p:p['replication']['pilot_seeds'].append(1000)),
               ('bool',lambda p:p['feeding'].update(contrast_nodes=True))]
        for name,mutate in cases:
            p=copy.deepcopy(self.protocol);mutate(p)
            with self.subTest(name=name),self.assertRaises(ValueError):validate_protocol(p)

    def test_json_roundtrip_and_restoration_preserve_every_state_field(self):
        u=arm_inputs(self.pair,'U')
        restored=json.loads(json.dumps(u))
        sim=build_simulation(restored,species_order=('resident','focal'))
        by_id={b.agent_id:b for sp in sim.species for b in sp.members}
        for record in restored['founders']:
            b=by_id[record['agent_id']]
            state=json.loads(json.dumps({k:v for k,v in b.__dict__.items() if k!='species'}))
            # The ledger may add no ecological state until the first tick.
            self.assertEqual(state,record['state'])
            self.assertEqual(b.species.population_role,record['population'])
        rng_before=random.getstate()
        restore_founders(sim,restored['founders'])
        self.assertEqual(random.getstate(),rng_before)
        self.assertTrue(sim._populated)

    def test_refined_clock_and_short_diagnostic_are_explicit(self):
        p=copy.deepcopy(self.protocol)
        p['clock'].update(dt_model_seconds=.05,horizon_model_seconds=2,steps=40)
        pair=make_pair_inputs(p,9)
        self.assertEqual(len(pair['feed_schedule']),2)
        self.assertEqual(set(tick_schedule(pair['feed_schedule'],.05)),{0,20})
        self.assertEqual(exact_steps(600,.05),12000)


if __name__=='__main__':unittest.main()
